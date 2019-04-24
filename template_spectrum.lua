package.path = package.path .. ";" .. arg[0]:gsub("[^/]*.lua","?.lua")
package.path = package.path .. ";./hmc/?.lua"
require 'common'
require 'smear'
require 'run'

-- Set number of colors to 4
qopqdp.defaultNc(4)
local Nc = 4
--profile = 1

local nf = 4
local nx = 24
local nt = 48
local mass = massGoesHere

-- Source and CG stuff
local src_start = 2       -- Time-slice of the first source
local src_num = 6         -- Nt/8 sources for Nt=48
local cg_prec = 1e-6      -- Max residual of CG
local cg_max = 5000       -- Maximum number of CG steps

-- Gauge fixing stuff
local gfix_prec = 1e-7    -- Max residual
local gfix_max = 5000     -- Max gauge fixing steps
local gfix_or = 1.75      -- Over-relaxation param for gauge fixing

-- Configuration to analyze
local inlat = 'inlatGoesHere'
printf("inlat = %s\n", inlat)

-- Use default mixed-precision inverter
local prec = 2

--- end of parameters

-- Print basic info and start preparing to load config
local latsize = { nx, nx, nx, nt }
local vol = 1
local spatvol = nx*nx*nx
local seed = seedGoesHere
printf("latsize =")
for k,v in ipairs(latsize) do
  vol=vol*v
  printf(" %i",v)
end
printf("\nvolume = %i\n", vol)
printf("mass = %g\n", mass)
printf("seed = %i\n", seed)
printf("gfix_prec = %g\n", gfix_prec)
printf("cg_prec = %g\n", cg_prec)
printf("src_num = %i\n", src_num)

-- Set up qopqdp
qopqdp.lattice(latsize);
qopqdp.profile(profile or 0);
qopqdp.verbosity(0);
qopqdp.seed(seed);

-- Start a timer
totaltime = qopqdp.dtime()

-- Load config
g = qopqdp.gauge()
if inlat then g:load(inlat)
else
  printf("No input configuration specified!\n")
  return 400
end

-- Reunitarize, just to be safe
do
  local devavg,devmax = g:checkSU()
  printf("unitarity deviation avg: %g  max: %g\n", devavg, devmax)
  g:makeSU()
  devavg,devmax = g:checkSU()
  printf("new unitarity deviation avg: %g  max: %g\n", devavg, devmax)
end

-- Print basic info about the config
function getplaq(g)
  local ps,pt = g:action{plaq=1}
  local lat = qopqdp.lattice()
  local nd,vol = #lat,1
  for i=1,nd do vol=vol*lat[i] end
  local s = 0.25*nd*(nd-1)*vol
  printf("plaq ss: %-8g  st: %-8g  tot: %-8g\n", ps/s, pt/s, 0.5*(ps+pt)/s)
end
getplaq(g)

-- Gauge fix to Coulomb gauge (0->x, ..., 3->t)
local t0 = qopqdp.dtime()
g:coulomb(3, gfix_prec, gfix_max, gfix_or)
t0 = qopqdp.dtime() - t0
printf("Coulgauge meas time: %g\n", t0)
getplaq(g)        -- Sanity check

-- Set up and run HYP smearing (verified to be consistent with MILC)
local smear = {}
smear[#smear+1] = { type="hyp", alpha={0.4,0.5,0.5} }
myprint("smear = ", smear, "\n")

-- We need to do this because 'smearGauge' expects an action object
printf("Start smearing\n")
local sg = smearGauge({g = g}, smear)
printf("Smearing done\n")
getplaq(sg)       -- Sanity check

-- Set the asqtad coeffs corresponding to plain staggered
coeffs = { one_link=1 }
w = qopqdp.asqtad();
w:coeffs(coeffs);
w:set(sg, prec);

-- Function to do an inversion and print out timing information
-- Based on actmt.set in asqtadact.lua
function solve_printinfo(w, dest, src, m, res, sub, opts)
  local t0 = qopqdp.dtime();
  w:solve({dest}, src, {m}, res, sub, opts);
  local cgtime = qopqdp.dtime() - t0;
  local flops = w:flops();
  local its = w:its();

  printf("inversion its: %g  secs: %g  Mflops: %g\n", its, cgtime, flops*1e-6/cgtime);
end

-- Function to add a correlator to an existing correlator array
-- Properly shift, normalize, and specialize to mesonic target
-- corr: What we're adding a new piece of data to
-- new_data: The new correlator we're adding in
-- t_shift: The 't' value of the source
-- norm: The normalization scaling new_data before adding it
function add_correlator_meson(corr, new_data, t_shift, norm)
  for j=1,#new_data do
    -- Compensate for shifted wall source
    j_real = (j+t_shift-1)%(#new_data)+1
    corr[j] = new_data[j_real]*norm + (corr[j] or 0)
  end
end

-- Function to build staggered phases
-- Use 0 for no phase, 1 for phase.
function make_phase_term(xsign,ysign,zsign,tsign)
  return xsign + 2*ysign + 4*zsign + 8*tsign
end

-- Set up where to put sources.
-- Evenly space src_num sources starting at t=src_start
local time_sources = {};
for i=1,src_num do
  time_sources[i] = (src_start+(i-1)*math.floor(nt/src_num))%nt;
end

-- Measure the pion using a random wall source
-- Equivalent to point source in the infinite stochastic source limit
-- We want Nc random sources
num_rand = Nc
norm = 1.0/(#time_sources*num_rand*spatvol*spatvol)

-- Need two "quark" objects
-- One holds the source, the other the sink
src = w:quark()
src:zero()
dest = w:quark()
dest:zero()

t = {}            -- A temporary place for correlator output
pion_fpi_ptp = {} -- Point to point fpi output
pion_fpi_ptw = {} -- Point to wall fpi output

do
  -- Loop over all time sources and all stochastic sources
  for srcnum=1,#time_sources do
    for randnum=1,num_rand do
      printf("Start random wall source %i at t=%i.\n", randnum, time_sources[srcnum]);
      io.stdout:flush();

      -- Create a random gaussian wall source
      src:wall_gaussian(time_sources[srcnum]);
      printf("gauss_src norm2: %g\n", src:norm2());

      -- And invert! This is a special function that also prints out CG info
      solve_printinfo(w, dest, src, mass, cg_prec, "all", {prec = prec, restart = cg_max});

      -- Contract for point sink
      t = dest:norm2("timeslices");
      add_correlator_meson(pion_fpi_ptp, t, time_sources[srcnum], norm);

      -- Contract for wall sink
      t = dest:norm2_wallsink("timeslices");
      add_correlator_meson(pion_fpi_ptw, t, time_sources[srcnum], norm);
    end
  end
end

-- Next, wall-source solves will need more output "quarks"
even_src = src; even_src:zero()
odd_src = dest; even_src:zero()
odd_soln = w:quark(); odd_soln:zero()
o_gupta = w:quark(); o_gupta:zero()
Do_gupta = w:quark(); Do_gupta:zero()
Dq_gupta = w:quark(); Dq_gupta:zero()
temp1 = w:quark(); temp1:zero()
temp2 = w:quark(); temp2:zero()
-- Nc-component arrays
even_soln, q_gupta = {}, {};
for i=1,Nc do
  even_soln[i] = w:quark(); even_soln[i]:zero();
  q_gupta[i] = w:quark(); q_gupta[i]:zero();
end

-- Since we have gauge fixed,
-- we don't need to use the gauge field in the parallel transporter
local unitg = qopqdp.gauge();
unitg:unit();

-- More destinations for results
p5, p5_g4, pion_ps_ck, pion_4_ck, pion_i5, pion_ij = {},{},{},{},{},{}
rho_0, rho_is, rho_ij, rho_i5 = {}, {}, {}, {}

-- Go!
norm = 1.0/(#time_sources)
do
  -- Loop over all time sources
  for srcnum=1,#time_sources do
    printf("Start wall source %i at t=%i.\n", srcnum, time_sources[srcnum]);

    -- Meson measurements, looping over Nc
    for i=1,Nc do
      printf("Start color %i.\n", i);
      io.stdout:flush();

      -- Prepare the even source: Even wall, norm matches MILC
      even_src:zero();
      even_src:wall(time_sources[srcnum], 0, i, -0.125);
      printf("even_src norm2: %g\n", even_src:norm2());

      -- Prepare the odd source: Odd wall, norm matches MILC
      odd_src:zero();
      odd_src:wall(time_sources[srcnum], 1, i, -0.125);
      printf("odd_src norm2: %g\n", odd_src:norm2());

      -- Invert on the even source
      -- We save each even solution even though we no longer compute the nucleon
      solve_printinfo(w, even_soln[i], even_src, mass, cg_prec, "all", {prec = prec, restart = cg_max});

      -- Invert on the odd source
      solve_printinfo(w, odd_soln, odd_src, mass, cg_prec, "all", {prec = prec, restart = cg_max});

      -- Follow MILC local pion computation
      -- These are the Goldstone pion, and the extra gamma_4 pion
      -- Only printed for a consistency check
      do -- Pion 0_A^(-+)
        t = even_soln[i]:norm2("timeslices");
        add_correlator_meson(pion_ps_ck, t, time_sources[srcnum], norm);
      end

      do -- Pion 0_A^(-+) w/ extra gamma_4
        temp1:set(even_soln[i]); -- Make a copy for rephasing
        temp1:rephase(make_phase_term(1,1,1,0), {0,0,0,time_sources[srcnum]}) -- (-1)^(x+y+z)
        t = even_soln[i]:Re_dot(temp1, "timeslices");
        add_correlator_meson(pion_4_ck, t, time_sources[srcnum], norm);
      end

      -- In the language of Gupta et al., the even_soln is "o+q"
      --                                  The  odd_soln is "q-o"
      -- We can thus reconstruct Gupta's o, q
      q_gupta[i]:zero();
      o_gupta:zero();
      q_gupta[i]:combine({even_soln[i], odd_soln},{1.0,1.0});
      o_gupta:combine({even_soln[i], odd_soln},{1.0,-1.0});

      -- Next prepare symmetric shifted values in the z direction
      -- We then zero out odd z coordinates
      -- These are in the wrong place with respect to the 2^4 hypercube,
      -- and will contribute to a wrong answer

      -- First, symm shift the odd solution in the z-direction
      -- Then remove odd z values
      Do_gupta:symshift(o_gupta, unitg, Nc);
      temp1:set(Do_gupta);
      temp1:rephase(make_phase_term(0,0,1,0),{0,0,0,time_sources[srcnum]})
      temp2:set(Do_gupta);
      Do_gupta:combine({temp1, temp2}, {0.5,0.5});

      -- Do the same to the even solution
      Dq_gupta:symshift(q_gupta[i], unitg, Nc);
      temp1:set(Dq_gupta);
      temp1:rephase(make_phase_term(0,0,1,0),{0,0,0,time_sources[srcnum]})
      temp2:set(Dq_gupta);
      Dq_gupta:combine({temp1, temp2}, {0.5,0.5});

      -- Using these states, we can compute the one-link-separated pion_i5 and pion_ij
      -- First, pion_i5 requires no phasing, just q Dq - o Do.
      do
        -- q Dq
        t = q_gupta[i]:Re_dot(Dq_gupta, "timeslices");
        add_correlator_meson(pion_i5, t, time_sources[srcnum], norm);

        -- minus o Do
        t = o_gupta:Re_dot(Do_gupta, "timeslices");
        add_correlator_meson(pion_i5, t, time_sources[srcnum], -norm);
      end

      -- Next, pion_ij requires rephasing by a factor of (-1)^(x+y+z), then o Dq - q Do
      Dq_gupta:rephase(make_phase_term(1,1,1,0), {0, 0, 0, time_sources[srcnum]});
      Do_gupta:rephase(make_phase_term(1,1,1,0), {0, 0, 0, time_sources[srcnum]});
      do
        -- q Dq
        t = o_gupta:Re_dot(Dq_gupta, "timeslices");
        add_correlator_meson(pion_ij, t, time_sources[srcnum], norm);

        -- minus o Do
        t = q_gupta[i]:Re_dot(Do_gupta, "timeslices");
        add_correlator_meson(pion_ij, t, time_sources[srcnum], -norm);
      end

      -- Better values for the pion and the g4_pion
      -- Taking combinations of q, o reduces errors
      -- The Goldstone pion is q q + o o
      do
        -- q q
        t = q_gupta[i]:Re_dot(q_gupta[i], "timeslices");
        add_correlator_meson(p5, t, time_sources[srcnum], norm);

        -- plus o o
        t = o_gupta:Re_dot(o_gupta, "timeslices");
        add_correlator_meson(p5, t, time_sources[srcnum], norm);
      end

      -- Next, the gamma4 pion is q o with a (-1)^(x+y+z) phase factor
      temp1:set(o_gupta);
      temp1:rephase(make_phase_term(1,1,1,0), {0, 0, 0, time_sources[srcnum]})
      do
        t = q_gupta[i]:Re_dot(temp1, "timeslices");
        add_correlator_meson(p5_g4, t, time_sources[srcnum], norm);
      end

      -- Shifting in the x, y directions gets us some rhos
      -- First, shifting in x gives
      --   1) rho_0: gamma_1 gamma_4 x taste_4
      --   2) rho_is: gamma_1 x 1

      -- Odd soln
      Do_gupta:symshift(o_gupta, unitg, 1);
      temp1:set(Do_gupta);
      temp1:rephase(make_phase_term(1,0,0,0),{0,0,0,time_sources[srcnum]})
      temp2:set(Do_gupta);
      Do_gupta:combine({temp1, temp2}, {0.5,0.5});

      -- Do the same to the even solution
      Dq_gupta:symshift(q_gupta[i], unitg, 1);
      temp1:set(Dq_gupta);
      temp1:rephase(make_phase_term(1,0,0,0),{0,0,0,time_sources[srcnum]})
      temp2:set(Dq_gupta);
      Dq_gupta:combine({temp1, temp2}, {0.5,0.5});

      -- First, rho_0 requires no phasing, just q Dq - o Do
      do
        -- q Dq
        t = q_gupta[i]:Re_dot(Dq_gupta, "timeslices");
        add_correlator_meson(rho_0, t, time_sources[srcnum], norm);

        -- minus o Do
        t = o_gupta:Re_dot(Do_gupta, "timeslices");
        add_correlator_meson(rho_0, t, time_sources[srcnum], -norm);
      end

      -- Next, rho_is requires rephasing by a factor of (-1)^(x+y+z), then o Dq - q Do
      Dq_gupta:rephase(make_phase_term(1,1,1,0), {0, 0, 0, time_sources[srcnum]});
      Do_gupta:rephase(make_phase_term(1,1,1,0), {0, 0, 0, time_sources[srcnum]});
      do
        -- q Dq
        t = o_gupta:Re_dot(Dq_gupta, "timeslices");
        add_correlator_meson(rho_is, t, time_sources[srcnum], norm);

        -- minus o Do
        t = q_gupta[i]:Re_dot(Do_gupta, "timeslices");
        add_correlator_meson(rho_is, t, time_sources[srcnum], -norm);
      end

      -- Next, shifting in y gives
      --   1) rho_ij: gamma_3 x taste_1 taste_4 taste_5
      --   2) rho_i5: gamma_3 gamma_4 x taste_1 taste_5

      -- Odd soln
      Do_gupta:symshift(o_gupta, unitg, 2);
      temp1:set(Do_gupta);
      temp1:rephase(make_phase_term(0,1,0,0),{0,0,0,time_sources[srcnum]})
      temp2:set(Do_gupta);
      Do_gupta:combine({temp1, temp2}, {0.5,0.5});

      -- Do the same to the even solution
      Dq_gupta:symshift(q_gupta[i], unitg, 2);
      temp1:set(Dq_gupta);
      temp1:rephase(make_phase_term(0,1,0,0),{0,0,0,time_sources[srcnum]})
      temp2:set(Dq_gupta);
      Dq_gupta:combine({temp1, temp2}, {0.5,0.5});

      -- First, rho_ij requires no phasing, just q Dq - o Do.
      do
        -- q Dq
        t = q_gupta[i]:Re_dot(Dq_gupta, "timeslices");
        add_correlator_meson(rho_ij, t, time_sources[srcnum], norm);

        -- minus o Do
        t = o_gupta:Re_dot(Do_gupta, "timeslices");
        add_correlator_meson(rho_ij, t, time_sources[srcnum], -norm);
      end

      -- Next, rho_is requires rephasing by a factor of (-1)^(x+y+z), then o Dq - q Do.
      Dq_gupta:rephase(make_phase_term(1,1,1,0), {0, 0, 0, time_sources[srcnum]});
      Do_gupta:rephase(make_phase_term(1,1,1,0), {0, 0, 0, time_sources[srcnum]});
      do
        -- q Dq
        t = o_gupta:Re_dot(Dq_gupta, "timeslices");
        add_correlator_meson(rho_i5, t, time_sources[srcnum], norm);

        -- minus o Do
        t = q_gupta[i]:Re_dot(Do_gupta, "timeslices");
        add_correlator_meson(rho_i5, t, time_sources[srcnum], -norm);
      end

      printf("End color\n", i);
    end

    printf("End wall source\n");
  end
end

-- Print results
printf("BEGIN_SPECTRUM\n");

--Reproduce MILC output format for POINT_KAON_5 WALL_KAON_5
printf("STARTPROP\n");
printf("MASSES:\t%.6e\t%.6e\n", mass, mass);
printf("SOURCE: RANDOM_WALL\n");
printf("SINKS: POINT_KAON_5 WALL_KAON_5\n");

for i = 1,#(pion_fpi_ptp) do
  printf("%i %.6e %f %.6e %f\n", i-1, pion_fpi_ptp[i], 0.0, pion_fpi_ptw[i], 0.0);
end
printf("ENDPROP\n");

--Reproduce MILC output format for PION_PS PION_SC
printf("STARTPROP\n");
printf("MASSES:\t%.6e\t%.6e\n", mass, mass);
printf("SOURCE: EVEN_WALL\n");
printf("SINKS: PION_PS PION_SC\n");

for i = 1,#(pion_ps_ck) do
  printf("%i %.6e %f %.6e %f\n", i-1, pion_ps_ck[i], 0.0, pion_4_ck[i], 0.0);
end
printf("ENDPROP\n");

--Reproduce MILC output format for PION_PS PION_SC PION_i5 PION_ij
printf("STARTPROP\n");
printf("MASSES:\t%.6e\t%.6e\n", mass, mass);
printf("SOURCE: EVENANDODD_WALL\n");
printf("SINKS: PION_PS PION_SC PION_i5 PION_ij\n");

for i = 1,#(p5) do
  printf("%i %.6e %f %.6e %f %.6e %f %.6e %f\n", i-1, p5[i], 0.0, p5_g4[i], 0.0, pion_i5[i], 0.0, pion_ij[i], 0.0);
end
printf("ENDPROP\n");

--Reproduce MILC output format for RHO_0 RHO_is RHO_ij RHO_i5
printf("STARTPROP\n");
printf("MASSES:\t%.6e\t%.6e\n", mass, mass);
printf("SOURCE: EVENANDODD_WALL\n");
printf("SINKS: RHO_0 RHO_is RHO_ij RHO_i5\n");

for i = 1,#(rho_0) do
  printf("%i %.6e %f %.6e %f %.6e %f %.6e %f\n", i-1, rho_0[i], 0.0, rho_is[i], 0.0, rho_ij[i], 0.0, rho_i5[i], 0.0);
end
printf("ENDPROP\n");

-- Done
printf("END_SPECTRUM\n");

totaltime = qopqdp.dtime() - totaltime;
printf("Total time: %f seconds.\n", totaltime);
io.stdout:flush()

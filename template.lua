package.path = package.path .. ";" .. arg[0]:gsub("[^/]*.lua","?.lua")
package.path = package.path .. ";./hmc/?.lua"
require 'common'
require 'run'

-- Set number of colors to 4
qopqdp.defaultNc(4)
--profile = 1

local nf = 4
local nx = 8
local nt = 4
local mass = massGoesHere
local beta = betaGoesHere
local beta_a = -0.25
local u0 = 1
local mass2 = mass2GoesHere
local prec = 2
local faresid = 1e-6
local grresid = 1e-6
local mdresid = 1e-6
local restart = 5000
--printf("faresid: %g\n", faresid)

local rhmc = {}
local hmcmasses = { mass, mass2 }

local outlat = nil
local ntraj = ntrajGoesHere
local tau = 1
local traj = trajGoesHere
local nsteps = nstepsGoesHere
local nsteps2 = nsteps2GoesHere
local ngsteps = ngstepsGoesHere
local nfsteps = { nsteps, nsteps2 }
nfsteps = repelem(nfsteps, nf/4)
local grcg = { prec=prec, resid=grresid, restart=restart }
local facg = { prec=prec, resid=faresid, restart=restart }
local mdcg = { prec=prec, resid=mdresid, restart=restart }
local ffprec = prec
local gintalg = {type="omelyan", lambda=lambdaG}
local fintalg = {type="omelyan", lambda=lambdaF}

-- Print out some useful data that don't seem to be recorded elsewhere
printf("ntraj %d\n", ntraj)
printf("traj_length %d\n", tau)
printf("Hasenbusch_mass %g\n", mass2)

local pbp = {}
pbp[1] = { reps=1 }
pbp[1].mass = mass
pbp[1].resid = 1e-6
pbp[1].opts = { restart=5000, max_restarts=1, max_iter=5000 }

-- Comment out next line to start from unit configuration
local latpat = 'latpatGoesHere'

--- end of parameters

function setpseudo(rhmc, mf, mb)
  local s1 = 4*mf*mf
  if mb then -- term (A+4mb^2)/(A+4mf^2)
    local s2 = 4*mb*mb
    local sr = math.sqrt(s1*s2)
    local c1 = s1 + s2 - 2*sr
    rhmc[#rhmc+1] = {
      GR = {{ {2, s2}, allfaceven=-2*mf, allfacodd=1, allmass2=-mb }},
      FA = {{ {1} }, { {math.sqrt(s2-s1), s1}, allfaceven=2*mf, allfacodd=1 }},
      MD = { {s2-s1, s1} }
    }
  else -- term 1/(A+4mf^2)
    rhmc[#rhmc+1] = {
      GR = {{ {1}, allfaceven=2*mf, allfacodd=-1 }},
      FA = {{ {1, s1}, allfaceven=2*mf, allfacodd=1 }},
      MD = { {1, s1} }
    }
  end
end
for i=1,#hmcmasses do
  for j=1,nf/4 do
    if i<#hmcmasses then
      setpseudo(rhmc, hmcmasses[i], hmcmasses[i+1])
    else
      setpseudo(rhmc, hmcmasses[i])
    end
  end
end
local npseudo = #rhmc
--printf("npseudo: %d\n", npseudo)

local p = {}
p.latsize = { nx, nx, nx, nt }
p.seed = seedGoesHere
p.beta = beta
p.nf = nf
p.u0 = u0
p.gaugeact = {type="plaquette_adjoint", adjFac=beta_a}
p.npseudo = npseudo

local smear = {}
smear[#smear+1] = { type="hyp", alpha={0.4,0.5,0.5} }
--myprint("smear = ", smear, "\n")
coeffs = { one_link=1 }
p.fermact = {type="asqtad", smear=smear, coeffs=coeffs, rhmc=rhmc}

local acts = setupacts(p)

local r = {}
r.ntraj = ntraj
r.tau = tau
local fp = {}
r.forceparams = fp
fp[1] = {}
fp[1][1] = {nsteps=ngsteps, intalg=gintalg}
local rhmc1 = {}
for j=1,npseudo do
  --printf("j = %i\n", j)
  --printf("ns = %i\n", nfsteps[j])
  fp[1][j+1] = {nsteps=nfsteps[j], intalg=fintalg}
  rhmc1[j] = {GR={},FA={},MD={}}
  rhmc1[j].GR[1] = {}
  rhmc1[j].GR[1].resid = grcg.resid
  rhmc1[j].FA.resid = facg.resid
  rhmc1[j].MD.resid = mdcg.resid
  rhmc1[j].GR[1].solveopts = {
    prec = grcg.prec,
    restart = grcg.restart
  }
  rhmc1[j].FA.solveopts = {
    prec = facg.prec,
    restart = facg.restart
  }
  rhmc1[j].MD.solveopts = {
    prec = mdcg.prec,
    restart = mdcg.restart,
    use_prev_soln = 1
  }
  rhmc1[j].MD.ffprec = ffprec
end
--myprint("rhmc1 = ", rhmc1, "\n")
copyto(rhmc, rhmc1)
--myprint("rhmc = ", rhmc, "\n")

r.pbp = pbp
--myprint("runparams = ", r, "\n")

if inlat then
  acts:load(inlat)
else
  if latpat then
    inlat = string.format(latpat, traj)
    acts:load(inlat)
  else
    acts:unit()
  end
end

local ps,pt = acts.fields.G:plaq()
printf("plaq ss: %g  st: %g  tot: %g\n", ps, pt, 0.5*(ps+pt))

r.md = false    --activate accept/reject step
r.ntraj = ntraj
acts:run(r)

-- Use plaq as stamp to check run
traj = traj + ntraj
local ps,pt = acts.fields.G:plaq()
printf("plaq ss: %g  st: %g  tot: %g\n", ps, pt, 0.5*(ps+pt))

-- Uncomment next line if starting from unit configuration above
--local latpat = 'latpatGoesHere'
if latpat then
  outlat = string.format(latpat, traj)
end
if outlat then
  acts:save(outlat)
end

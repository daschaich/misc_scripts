package.path = package.path .. ";" .. arg[0]:gsub("[^/]*.lua","?.lua")
package.path = package.path .. ";./hmc/?.lua"
require 'common'
require 'run'

--profile = 1

--u0 is needed by setupacts()
local nx = 48
local nt = 96
local u0 = 1
local prec = 2

local rhmc = {}

local p = {}
p.latsize = { nx, nx, nx, nt }
p.seed = 41
p.beta = 4.8
p.nf = 8
p.u0 = u0
p.gaugeact = {type="plaquette_adjoint", adjFac=-0.25}
p.npseudo = 1
coeffs = { one_link=1 }
p.fermact = {type="asqtad", coeffs=coeffs, rhmc=rhmc}

local acts = setupacts(p)
local inlat = 'inlatGoesHere'

local outlat = 'outlatGoesHere'
acts:load(inlat)
acts:save(outlat)

using Random
using ParameterizedFunctions
using DifferentialEquations
using ODE
using ODEInterfaceDiffEq
using Plots
using JLD
using GLMNet
using Latexify

# loops over the tesselation model using different random random initial conditions each time,
# then, saves the output as a matlab mat file in the home directory of the simulation with the following name
# <sim dir>/[<original summary mat file name>,<julia_sim>].mat

# for each species, the data is saved in arrays: u and du, each is a |time| x |cells| matrix
root_dir = "/home/au/mkarikom/Software/julia_voronoi_iSINDy/"
#root_dir = "/home/au/Insync/mkarikom@uci.edu/Google Drive/Software/julia_voronoi_iSINDy/"

julia_dir = "./Julia_scripts"
tesselation = "Simulation_Voronoi.17_T.1000_Samplerate.1__09-27-20-16-16"
tesselation_file = string("./Sim_Summary/", tesselation,".mat")
result_file = string("./Sim_Summary/",tesselation, "_juliaSim_Linear_Loop")

include("runLattice.jl")
include("runLatticeSmall.jl")

cd(root_dir)

using Pkg
Pkg.activate(julia_dir)

using MAT
# read in the data for the adjacency and tesselation
#file = matread("/home/au/Insync/mkarikom@uci.edu/Google Drive/Software/julia_voronoi_iSINDy/Sim_Summary/Simulation_Voronoi.17_T.1000_Samplerate.1__09-26-13-01-34.mat")
file = matread(tesselation_file)

dist = file["distances"]
adj = file["adjacency"]
n_vars = 2 # this is the number of molecular species e.g. |{notch, delta}| = 2 for the simple model
n_runs = 10
tspan = (0.0,30.0)
mult = .01 # the coefficient for the gaussian measurement noise
trim_ends = 15 # trim the ends to get rid of weird interpolation issues
Random.seed!(123)
n_cells = size(dist)[1]
myalpha = 1


p = Array{Float64}(undef, n_cells*5)
for pp = 1:n_cells
    p[(pp-1)*5+1] = 1 # nu
    p[(pp-1)*5+2] = 5 # betaD
    p[(pp-1)*5+3] = 20 # betaR
    p[(pp-1)*5+4] = 3 # h
    p[(pp-1)*5+5] = 3 # m
end

perf = DataFrame()
check = runPerturb(p, n_cells, trim_ends, mult, myalpha, tspan, n_runs, n_vars, adj, dist)

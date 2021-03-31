# Test case with simple looping flow field
# streamfunction is given by s=sin(pi*x)*sin(pi*y)

using Particles
using Plots
using Random

randpool = MersenneTwister(0)

# collected configuration is in Dict d
d = default_userdata() # start with some defaults
# settings for this experiment
n = 30 # number of particles
d["nparticles"] = n
# all variables for one particle are collected in a vector
variables = ["x", "y", "z", "age"]
d["variables"] = variables
# initial position of the particles
m = length(variables)
p = zeros(m, n)
p[1, :] = 0.2 .+ 0.02 * randn(randpool, n, 1)
p[2, :] = 0.2 .+ 0.02 * randn(randpool, n, 1)
d["particles"] = p # initial values
# simulation time
d["dt"] = 0.01
d["tstart"] = 0.0
d["tend"] = 2.5
# write to netcdf
d["write_maps_times"] = collect(0.0:0.1:2.5)
d["write_maps"] = false
d["write_maps_filename"] = "output_loop.nc"
# write plots to file
d["plot_maps_times"] = collect(0.0:0.5:2.5)
d["plot_maps"] = false

# forcing currents stream (used only for plotting here) , u and v
function stream(x, y, z, t)
    sin(pi * x) * sin(pi * y)
end
# flow in x direction (for now has to be called u)
function u(x, y, z, t)
   # u=-s_y
    -sin(pi * x) * cos(pi * y)
end

# flow in y direction (for now has to be called v)
function v(x, y, z, t)
   # v=s_x
    cos(pi * x) * sin(pi * y)
end


"""
   !f(ds,s,t,i,d)

Dynamic model, computes time derivative ds of s at current time t
for particle i and possibly using data/functions from d of type userdata.
"""
function f!(∂s, s, t, i, d)
   x, y, z, age = s
   # dx/dt=u
    ∂s.x = u(x, y, z, t)
   # dy/dt=v
    ∂s.y = v(x, y, z, t)
   # dz/dt=0
    ∂s.z = 0.0
   # age=(t-t0)
    ∂s.t = 1.0
end
d["f"] = f!

function plot_background(d)
   # streamfunction for plot
    x1 = 0.0:0.01:1.0
    y1 = 0.0:0.01:1.0
    flow = zeros(length(x1), length(y1))
    for i = 1:length(x1)
        for j = 1:length(y1)
            flow[i,j] = stream(x1[i], y1[j], 0.0, 0.0)
        end
    end
    f = contour(x1, y1, flow', legend=false)
    return(f)
end
d["plot_maps_background"] = plot_background

println("Start with run_simulation(d) if it does not start automatically")
# run_simulation(d)  # compile
# @time run_simulation(d)
nothing

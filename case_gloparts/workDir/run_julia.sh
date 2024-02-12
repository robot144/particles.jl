#! /bin/bash
export PATH=$PATH:/opt/fews/models/julia-1.7.3/bin
julia --project= ../Particles.jl/drifters.jl ../workDir/config.toml > out.txt

# package Particles
Particle modelling tools for coastal flow

# using the package
- add your julia development folder to the LOAD_PATH. This can be done on the prompt or in your startup.jl file
push!(LOAD_PATH, "$(homedir())/src/julia")  #adapt path to your local folder structure $HOME/.julia/config/startup.jl
or push!(LOAD_PATH, "d:\\some_folder\\julia" #windows c:\Users\username\.julia\startup.jl
- add the package Particles
using Pkg
Pkg.add("Particles")
or activate the julia package management environment with the ']'-key on the REPL and use commands like
add Particles
# testing
Pkg.test("Particles")
# to use the package
using Particles

@ECHO OFF
ECHO This is to execute Julia Global Model in GLOSSIS

julia --project= ../Particles.jl/drifters.jl ../workDir/config.toml > out.txt

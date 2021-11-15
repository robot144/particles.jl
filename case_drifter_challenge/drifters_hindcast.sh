#! /bin/bash

## csv file allows simple control
#julia --project=.. drifters.jl config_hindcast.csv
## more complex input in toml files (https://toml.io/en/)
julia --project=.. drifters.jl config_hindcast.toml

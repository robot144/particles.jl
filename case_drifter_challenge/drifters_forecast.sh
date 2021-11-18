#! /bin/bash

## more complex input in toml files (https://toml.io/en/)
julia --project=.. drifters.jl config_forecast.toml

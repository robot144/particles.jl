using TOML
using CSV

function run_simulation()
    usage = "Usage: julia -e 'using Particles; run_simulation()' '/path/to/config.toml/csv'"
    n = length(ARGS)
    if n != 1
        throw(ArgumentError(usage))
    end
    config_path = only(ARGS)
    if !isfile(config_path)
        throw(ArgumentError("File not found: $(config_path)\n" * usage))
    end
    run_simulation(config_path)
end


function run_simulation(path::AbstractString)
    d = config(path)
    run_simulation(d)
end

function config(path::AbstractString)
    _, ext = splitext(path)
    if ext == ".csv"
        c = CSV.File(path; normalizenames=true)
        length(c) < 1 && error("Empty csv provided.")
        config = Dict(string.(keys(c[1])) .=> values(c[1]))
    elseif ext == ".toml"
        config = TOML.parsefile(path)
    else
        error("Unknown file format $ext provided.")
    end
    merge(default_userdata(), config)
end

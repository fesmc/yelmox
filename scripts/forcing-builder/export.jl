# export.jl — write docs/forcing-builder/forcing.json, the data of the forcing
# builder page (docs/forcing-builder.md): the schema (input/yelmox_defaults.nml
# with its comments), the forcing library (input/forcing_db.nml), the run
# control and forcing groups of the base par files (yelmox/yelmox_*.nml), the
# climate database groups (input/clim_db.nml, for display) and the scalar
# series files of input/ (downsampled), so the page can draw every tsgen
# series and the sea level, and the cases of the docs figures as presets
# (scripts/docs-figures/run.sh --list).
#
# Usage (from the repository root):
#   julia --project=scripts/forcing-builder scripts/forcing-builder/export.jl [OUT.json]

include(joinpath(@__DIR__, "Namelist.jl"))
include(joinpath(@__DIR__, "TsGen.jl"))
include(joinpath(@__DIR__, "Resolve.jl"))

using .Namelist, .TsGen, .Resolve
using JSON, Dates, NCDatasets

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = length(ARGS) >= 1 ? ARGS[1] : joinpath(ROOT, "docs", "forcing-builder", "forcing.json")
const NMAX = 800          # points per exported series

# The groups of a base exported: run control, components, forcing, transforms.
const BASE_GROUPS = r"^(ctrl|comps|sim|output|tm_\w+|forcing|coef_\w+|term_\w+|snapshots|snapshot_\w+|taper|bias|orbit|ice_bnd|barysealevel|esm|obm_forcing|domain|spinup|transient)$"

jkey(k::Key) = Dict("raw" => k.raw, "v" => parse_value(k.raw), "c" => k.comment, "pre" => k.pre)
jgroup(g::Group) = Dict("c" => g.comment, "order" => [k.name for k in g.keys],
                        "keys" => Dict(k.name => jkey(k) for k in g.keys))

function downsample(t::Vector{Float64}, f::Vector{Float64}, nmax::Int)
    n = length(t)
    n <= nmax && return (t, f)
    stride = ceil(Int, n / nmax)
    idx = collect(1:stride:n)
    idx[end] == n || push!(idx, n)
    return (t[idx], f[idx])
end

clean(x) = isfinite(x) ? x : nothing

function export_ascii(path)
    ser = load_series(path)
    t, f = downsample(ser.time, ser.f, NMAX)
    Dict("kind" => "ascii", "header" => ser.header, "nc" => ser.nc, "n" => length(ser.time),
         "time" => t, "f" => clean.(f))
end

function export_netcdf(path)
    vars = Dict{String,Any}()
    ds = NCDatasets.Dataset(path)
    try
        haskey(ds, "time") || return Dict("kind" => "netcdf", "vars" => vars)
        t = Float64.(collect(ds["time"][:]))
        for name in keys(ds)
            v = ds[name]
            name == "time" && continue
            NCDatasets.dimnames(v) == ("time",) || continue
            raw = collect(v[:])
            f = [x === missing ? NaN : Float64(x) for x in raw]
            mv = get(v.attrib, "missing_value", nothing)
            mv !== nothing && (f = [x == mv ? NaN : x for x in f])
            td, fd = downsample(t, f, NMAX)
            vars[name] = Dict("units" => string(get(v.attrib, "units", "")), "n" => length(t),
                              "time" => td, "f" => clean.(fd))
        end
    finally
        close(ds)
    end
    return Dict("kind" => "netcdf", "vars" => vars)
end

function main()
    cfg = load_config(ROOT)

    defaults = Dict(name => jgroup(g) for (name, g) in cfg.defaults)
    library = Dict("files" => cfg.library_files,
                   "groups" => Dict(name => jgroup(g) for (name, g) in cfg.library))

    bases = Dict{String,Any}()
    for f in sort(filter(x -> occursin(r"^yelmox_\w+\.nml$", x), readdir(joinpath(ROOT, "yelmox"))))
        name = replace(f, ".nml" => "")
        groups = read_namelist(joinpath(ROOT, "yelmox", f))
        sel = [g for g in groups if occursin(BASE_GROUPS, g.name)]
        bases[name] = Dict("file" => "yelmox/" * f, "order" => [g.name for g in sel],
                           "groups" => Dict(g.name => jgroup(g) for g in sel))
    end

    climdb = Dict{String,Any}()
    for g in read_namelist(joinpath(ROOT, "input", "clim_db.nml"))
        d = Dict(k.name => parse_value(k.raw) for k in g.keys if k.name in ("filename", "name", "time", "units", "grid_src"))
        climdb[g.name] = d
    end

    # series: every input/*.dat, plus the netCDF files any series_file/sl_path names
    series = Dict{String,Any}()
    for f in sort(readdir(joinpath(ROOT, "input")))
        endswith(f, ".dat") || continue
        p = "input/" * f
        try
            series[p] = export_ascii(joinpath(ROOT, p))
        catch err
            @warn "series $p skipped" err
        end
    end
    ncfiles = Set{String}()
    for src in (values(cfg.library), (read_namelist(joinpath(ROOT, b["file"])) for b in values(bases))...)
        for g in src, k in g.keys
            k.name in ("series_file", "sl_path") || continue
            v = parse_value(k.raw)
            v isa String && endswith(v, ".nc") && push!(ncfiles, v)
        end
    end
    for p in ncfiles
        isfile(joinpath(ROOT, p)) || continue
        try
            series[p] = export_netcdf(joinpath(ROOT, p))
        catch err
            @warn "series $p skipped" err
        end
    end

    # presets: the cases of the docs figures (scripts/docs-figures/run.sh --list),
    # each a base and the -p overrides of the experiment and its timeline
    presets = Any[]
    for line in eachline(`bash $(joinpath(ROOT, "scripts", "docs-figures", "run.sh")) --list`)
        cols = split(line, '\t')
        length(cols) >= 2 || continue
        base = replace(basename(cols[2]), ".nml" => "")
        haskey(bases, base) || (@warn "preset $(cols[1]): base $base not exported"; continue)
        push!(presets, Dict("name" => cols[1], "base" => base, "p" => String.(cols[3:end])))
    end

    commit = try
        strip(read(`git -C $ROOT rev-parse --short HEAD`, String))
    catch
        ""
    end
    data = Dict(
        "generated" => string(Dates.today()),
        "commit" => commit,
        "templates" => [Dict("pattern" => re.pattern, "template" => t) for (re, t) in Resolve.TEMPLATES],
        "defaults" => defaults,
        "library" => library,
        "bases" => bases,
        "climdb" => climdb,
        "series" => series,
        "presets" => presets,
    )
    mkpath(dirname(OUT))
    open(OUT, "w") do io
        JSON.print(io, data)
    end
    println("wrote $OUT ($(round(filesize(OUT)/1024)) kB): $(length(bases)) bases, $(length(cfg.library)) library groups, $(length(series)) series, $(length(presets)) presets")
end

main()

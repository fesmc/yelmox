# Export the Antarctica multigrid ensemble to a single compact JLD2 bundle.
#
# Reads the runs staged by yelmox/run_mg_resolution.sh and writes, per run, the
# final-time 2D maps + the full 1D timeseries, plus the PD observational targets
# (BedMachine topography, Rignot velocities) on each run's Yelmo grid. The result
# is a few MB (vs. hundreds of MB of NetCDF) so it can be copied to a laptop and
# analyzed/plotted locally with only JLD2 + CairoMakie (no NetCDF C library, no
# access to the raw output).
#
# By default only COMPLETED runs (reached ~25 kyr) are exported -- pass tmin<=0
# to export every run that has any output.
#
# Usage (from the yelmox root, on the cluster):
#   julia --project=analysis analysis/export_mg_jld2.jl [run_root] [out.jld2] [tmin]
#     run_root  default output/mg_opt
#     out.jld2  default analysis/mg_ensemble.jld2
#     tmin      default 24999   (min last-time [yr] to count as complete; <=0 = all)
#
# Load locally:
#   using JLD2
#   ens = load("mg_ensemble.jld2", "ens")
#   r   = ens["runs"]["y32KM_m32KM"]; r["maps"]["H_ice"]; r["ts"]["V_sle"]
#   o   = ens["obs"][r["ygrid"]];     o["H_ice"], o["uxy_s"]

using NCDatasets, JLD2, Printf

struct Run; name::String; ygrid::String; mgrid::String; end
runs = [
    Run("y8KM_m8KM",   "ANT-8KM",  "ANT-8KM"),
    Run("y16KM_m16KM", "ANT-16KM", "ANT-16KM"),
    Run("y32KM_m32KM", "ANT-32KM", "ANT-32KM"),
    Run("y32KM_m8KM",  "ANT-32KM", "ANT-8KM"),
    Run("y32KM_m16KM", "ANT-32KM", "ANT-16KM"),
    Run("y16KM_m8KM",  "ANT-16KM", "ANT-8KM"),
]

# Final-time 2D fields worth carrying (extracted only if present in the file).
MAP_VARS = ["H_ice", "z_srf", "z_bed", "uxy_s", "uxy_b", "uxy_bar",
            "bmb_shlf", "bmb", "smb", "mb_net", "cb_ref", "beta", "taub"]
# Full 1D timeseries worth carrying.
TS_VARS  = ["time", "V_ice", "V_sle", "V_ice_g", "V_ice_f", "A_ice", "A_ice_g",
            "A_ice_f", "z_srf", "uxy_s", "uxy_s_g", "bmb_shlf", "smb", "dVidt"]

yelmo2d(root, r) = joinpath(root, r.name, "yelmo_$(r.ygrid).nc")
yelmots(root, r) = joinpath(root, r.name, "yelmo_$(r.ygrid)_ts.nc")

# Locate the ice_data tree. Each run dir carries an `ice_data` symlink; fall back
# to a cwd / repo-root copy. Returns "" if none found.
function find_icedata(root, names)
    cands = String[]
    for n in names; push!(cands, joinpath(root, n, "ice_data")); end
    push!(cands, "ice_data", joinpath(root, "..", "..", "ice_data"))
    for c in cands; isdir(c) && return c; end
    return ""
end
obs_topo(base, g) = joinpath(base, "Antarctica", g, "$(g)_TOPO-BedMachine.nc")
obs_vel(base, g)  = joinpath(base, "Antarctica", g, "$(g)_VEL-R11-2.nc")

"Last time value in a 2D file, or -Inf if unreadable."
function run_tmax(root, r)
    p = yelmo2d(root, r)
    isfile(p) || return -Inf
    NCDataset(p) do ds
        (haskey(ds, "time") && length(ds["time"]) > 0) ? Float64(ds["time"][end]) : -Inf
    end
end

"Read final-time maps + grid from a run's 2D file."
function read_run_maps(path)
    NCDataset(path) do ds
        d = Dict{String,Any}()
        d["xc"] = Float64.(ds["xc"][:])
        d["yc"] = Float64.(ds["yc"][:])
        m = Dict{String,Array{Float32,2}}()
        for v in MAP_VARS
            haskey(ds, v) || continue
            m[v] = Float32.(ds[v][:, :, end])
        end
        d["maps"] = m
        return d
    end
end

"Read the full timeseries from a run's _ts file (empty Dict if absent)."
function read_run_ts(path)
    isfile(path) || return Dict{String,Vector{Float32}}()
    NCDataset(path) do ds
        t = Dict{String,Vector{Float32}}()
        for v in TS_VARS
            haskey(ds, v) || continue
            t[v] = Float32.(ds[v][:])
        end
        return t
    end
end

"Read PD obs (topography + velocity) on grid g; nothing if files missing."
function read_obs(base, g)
    tp, vp = obs_topo(base, g), obs_vel(base, g)
    (isfile(tp) && isfile(vp)) || return nothing
    o = Dict{String,Any}()
    NCDataset(tp) do ds
        o["xc"] = Float64.(ds["xc"][:]); o["yc"] = Float64.(ds["yc"][:])
        for v in ("H_ice", "z_bed", "z_srf")
            haskey(ds, v) && (o[v] = Float32.(ds[v][:, :]))
        end
    end
    NCDataset(vp) do ds
        v = haskey(ds, "uxy_srf") ? "uxy_srf" : (haskey(ds, "uxy_s") ? "uxy_s" : "")
        v != "" && (o["uxy_s"] = Float32.(ds[v][:, :]))
    end
    return o
end

function main()
    root = length(ARGS) >= 1 ? ARGS[1] : "output/mg_opt"
    out  = length(ARGS) >= 2 ? ARGS[2] : "analysis/mg_ensemble.jld2"
    tmin = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 24999.0

    runsout = Dict{String,Any}()
    grids   = String[]
    for r in runs
        tmax = run_tmax(root, r)
        if tmax < tmin
            @info @sprintf("skip %-14s (tmax=%.0f < tmin=%.0f)", r.name, tmax, tmin)
            continue
        end
        rd = read_run_maps(yelmo2d(root, r))
        rd["name"]  = r.name; rd["ygrid"] = r.ygrid; rd["mgrid"] = r.mgrid
        rd["tmax"]  = tmax
        rd["ts"]    = read_run_ts(yelmots(root, r))
        runsout[r.name] = rd
        r.ygrid in grids || push!(grids, r.ygrid)
        @info @sprintf("export %-14s  ygrid=%-8s mgrid=%-8s tmax=%.0f  maps=%d ts=%d",
                       r.name, r.ygrid, r.mgrid, tmax, length(rd["maps"]), length(rd["ts"]))
    end

    base = find_icedata(root, collect(keys(runsout)))
    base == "" && @warn "ice_data tree not found; obs will be empty"
    obs = Dict{String,Any}()
    for g in grids
        o = base == "" ? nothing : read_obs(base, g)
        o === nothing ? (@warn "obs missing for $g") : (obs[g] = o)
    end

    isempty(runsout) && (@warn "no runs met tmin=$tmin; nothing written"; return)

    ens = Dict("runs" => runsout, "obs" => obs,
               "meta" => Dict("run_root" => root, "tmin" => tmin,
                              "map_vars" => MAP_VARS, "ts_vars" => TS_VARS))
    mkpath(dirname(out) == "" ? "." : dirname(out))
    jldsave(out; ens = ens)
    sz = stat(out).size / 1e6
    @info @sprintf("wrote %s  (%d runs, %d obs grids, %.1f MB)",
                   out, length(runsout), length(obs), sz)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

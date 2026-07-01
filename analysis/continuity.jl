# Diagnostics for yelmox_mg runs:
#   (1) restart continuity  -- overlay a straight-through run and a run restarted
#       partway, per module, to show state is continuous across the restart seam.
#   (2) cross-driver parity -- overlay yelmox_mg against the reference yelmox.f90.
#
# Reads the climber-x per-module timeseries files (<module>_<grid>_ts.nc) written
# by yelmox_mg, and the yelmo1D.nc written by yelmox.f90.
#
# Usage (from the yelmox root):
#   julia --project=analysis analysis/continuity.jl <run_root> [out_dir]
#
# <run_root> is expected to contain:
#   straight/   a 0->T run (restart written partway)
#   restart/    a run restarted from the partway bundle
#   parity/ref/        yelmox.f90         (yelmo1D.nc)
#   parity/mg_id/      yelmox_mg identity (all grids = Yelmo grid)
#   parity/mg_hub16/   yelmox_mg hi-res hub
# Missing sub-runs are skipped.

using CairoMakie, NCDatasets, Printf

# --- reading helpers -------------------------------------------------------

"Find the module timeseries file <module>_*_ts.nc in dir (or a fallback name)."
function find_ts(dir, prefix; fallback=nothing)
    isdir(dir) || return nothing
    for f in readdir(dir)
        (startswith(f, prefix) && endswith(f, "_ts.nc")) && return joinpath(dir, f)
    end
    fallback !== nothing && isfile(joinpath(dir, fallback)) && return joinpath(dir, fallback)
    return nothing
end

"Read (time, var) from a netCDF file, or nothing if file/var absent."
function readts(path, var)
    (path === nothing || !isfile(path)) && return nothing
    NCDataset(path) do ds
        haskey(ds, var) || return nothing
        return (Float64.(ds["time"][:]), Float64.(ds[var][:]))
    end
end

yelmo_ts(dir) = find_ts(dir, "yelmo"; fallback="yelmo1D.nc")
isos_ts(dir)  = find_ts(dir, "isos")

# --- panels ----------------------------------------------------------------

"Overlay one variable from several runs onto an axis; mark the restart seam."
function panel!(ax, series; seam=nothing)
    for (label, path, var, kw) in series
        d = readts(path, var)
        d === nothing && continue
        lines!(ax, d[1], d[2]; label=label, kw...)
    end
    seam !== nothing && vlines!(ax, [seam]; color=:gray, linestyle=:dash)
end

# --- figures ---------------------------------------------------------------

function fig_continuity(root, out; seam=50.0)
    sdir = joinpath(root, "straight"); rdir = joinpath(root, "restart")
    sy, ry = yelmo_ts(sdir), yelmo_ts(rdir)
    si, ri = isos_ts(sdir),  isos_ts(rdir)

    specs = [
        ("Yelmo V_ice [1e6 km^3]", [("straight", sy, "V_ice", (;color=:black)),
                                    ("restart",  ry, "V_ice", (;color=:red, linestyle=:dash, linewidth=2))]),
        ("Yelmo A_ice [1e6 km^2]", [("straight", sy, "A_ice", (;color=:black)),
                                    ("restart",  ry, "A_ice", (;color=:red, linestyle=:dash, linewidth=2))]),
        ("isos bsl [m]",           [("straight", si, "bsl", (;color=:black)),
                                    ("restart",  ri, "bsl", (;color=:red, linestyle=:dash, linewidth=2))]),
        ("isos mean z_bed [m]",    [("straight", si, "z_bed_mean", (;color=:black)),
                                    ("restart",  ri, "z_bed_mean", (;color=:red, linestyle=:dash, linewidth=2))]),
        ("isos mean w [m]",        [("straight", si, "w_mean", (;color=:black)),
                                    ("restart",  ri, "w_mean", (;color=:red, linestyle=:dash, linewidth=2))]),
        ("isos mean we [m]",       [("straight", si, "we_mean", (;color=:black)),
                                    ("restart",  ri, "we_mean", (;color=:red, linestyle=:dash, linewidth=2))]),
    ]

    fig = Figure(size=(1000, 640))
    Label(fig[0, 1:3], "Restart continuity (restart at $(seam) yr)", fontsize=18)
    for (k, (ttl, series)) in enumerate(specs)
        r, c = fldmod1(k, 3)
        ax = Axis(fig[r, c]; title=ttl, xlabel="time [yr]")
        panel!(ax, series; seam=seam)
        k == 1 && axislegend(ax; position=:lb, framevisible=false)
    end
    save(out, fig)
    println("wrote ", out)
end

function fig_parity(root, out)
    pdir = joinpath(root, "parity")
    ref, mid, mhub = joinpath(pdir, "ref"), joinpath(pdir, "mg_id"), joinpath(pdir, "mg_hub16")
    specs = [
        ("Yelmo V_ice [1e6 km^3]", "V_ice"),
        ("Yelmo A_ice [1e6 km^2]", "A_ice"),
    ]
    fig = Figure(size=(900, 380))
    Label(fig[0, 1:2], "Parity: yelmox_mg vs yelmox.f90", fontsize=18)
    for (k, (ttl, var)) in enumerate(specs)
        ax = Axis(fig[1, k]; title=ttl, xlabel="time [yr]")
        panel!(ax, [
            ("yelmox.f90",   yelmo_ts(ref),  var, (;color=:black, linewidth=2)),
            ("mg identity",  yelmo_ts(mid),  var, (;color=:dodgerblue, linestyle=:dash, linewidth=2)),
            ("mg hub 16km",  yelmo_ts(mhub), var, (;color=:orange, linestyle=:dot, linewidth=2)),
        ])
        k == 1 && axislegend(ax; position=:lb, framevisible=false)
    end
    save(out, fig)
    println("wrote ", out)
end

# --- main ------------------------------------------------------------------

function main()
    root = length(ARGS) >= 1 ? ARGS[1] : error("usage: continuity.jl <run_root> [out_dir]")
    out  = length(ARGS) >= 2 ? ARGS[2] : joinpath(root, "figures")
    mkpath(out)
    isdir(joinpath(root, "straight")) && fig_continuity(root, joinpath(out, "restart_continuity.png"))
    isdir(joinpath(root, "parity"))   && fig_parity(root, joinpath(out, "parity.png"))
end

main()

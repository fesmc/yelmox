# Forcing figure of one climate-only run (run.sh): the climate at one or more
# forcing times minus the climate at a reference time (annual and summer air
# temperature difference, precipitation ratio, and the surface mass balance
# difference when the climate supplies it), or with --abs the climate at the
# reference time itself. The coast (z_srf = 0 of the climate file, else of
# yelmo.nc in the run folder: esm.nc has no z_srf) is drawn.
#
# Usage:
#   julia --project=scripts/docs-figures scripts/docs-figures/plot.jl RUN OUT.png \
#         --ref T [--at T1 T2 ...] [--ref-run RUN0] [--abs] [--title TEXT]
#
#   RUN    run folder (snap.nc, esm.nc or rembo.nc in it)
#   --ref  reference time (as in the file's time axis)
#   --at   forcing times to show minus the reference
#   --abs  show the reference climate itself (no --at)
#   --ref-run RUN0  take the reference from another run (e.g. the same par file
#          with the coefficient at present day)

using NCDatasets
using CairoMakie
using Printf

function parse_args(args)
    opt = Dict{String,Any}("run" => args[1], "out" => args[2], "ref" => NaN, "at" => Float64[],
                           "abs" => false, "title" => "", "ref_run" => "")
    i = 3
    while i <= length(args)
        a = args[i]
        if a == "--ref"
            opt["ref"] = parse(Float64, args[i+1]); i += 2
        elseif a == "--at"
            i += 1
            while i <= length(args) && !startswith(args[i], "--")
                push!(opt["at"], parse(Float64, args[i])); i += 1
            end
        elseif a == "--ref-run"
            opt["ref_run"] = args[i+1]; i += 2
        elseif a == "--abs"
            opt["abs"] = true; i += 1
        elseif a == "--title"
            opt["title"] = args[i+1]; i += 2
        else
            error("unknown option $a")
        end
    end
    return opt
end

function climate_file(run)
    for f in ("snap.nc", "esm.nc", "rembo.nc")
        p = joinpath(run, f)
        isfile(p) && return p
    end
    error("no climate file in $run")
end

# The coast: z_srf of the climate file at record k, else of the run's yelmo.nc
# (first record, on its own grid). Returns (x, y, z_srf) or nothing.
function coast(run, ds, k)
    haskey(ds, "z_srf") && return (Float64.(ds["xc"][:]), Float64.(ds["yc"][:]), field(ds, "z_srf", k))
    p = joinpath(run, "yelmo.nc")
    isfile(p) || return nothing
    NCDataset(p) do dy
        haskey(dy, "z_srf") ? (Float64.(dy["xc"][:]), Float64.(dy["yc"][:]), field(dy, "z_srf", 1)) : nothing
    end
end

field(ds, v, k) = haskey(ds, v) ? Float64.(coalesce.(ds[v][:, :, k], NaN)) : nothing

function record(t, t0)
    k = argmin(abs.(t .- t0))
    abs(t[k] - t0) > 1e-3 * max(1.0, abs(t0)) && error("time $t0 not in the file (nearest $(t[k]))")
    return k
end

tlabel(t) = abs(t) >= 1000 && t <= 0 ? @sprintf("%g ka", -t / 1000) : @sprintf("%g", t)

function main(opt)
    f  = climate_file(opt["run"])
    ds = NCDataset(f)
    x = Float64.(ds["xc"][:]); y = Float64.(ds["yc"][:]); t = Float64.(ds["time"][:])
    # the reference: this run, or another (--ref-run)
    ds0 = isempty(opt["ref_run"]) ? ds : NCDataset(climate_file(opt["ref_run"]))
    k0 = record(Float64.(ds0["time"][:]), opt["ref"])
    cst = coast(isempty(opt["ref_run"]) ? opt["run"] : opt["ref_run"], ds0, k0)
    aspect = (x[end] - x[1]) / (y[end] - y[1])

    if opt["abs"]
        panels = [("Air temperature, annual [°C]", field(ds0, "t2m_ann", k0) .- 273.15, (-40, 10), :thermal, identity),
                  ("Air temperature, summer [°C]", field(ds0, "t2m_sum", k0) .- 273.15, (-20, 20), :thermal, identity),
                  ("Precipitation [mm/a]", field(ds0, "pr_ann", k0), (50, 3000), :Blues, log10)]
        filter!(p -> p[2] !== nothing, panels)
        rows = [("", panels)]
    else
        rows = []
        for ta in opt["at"]
            k = record(t, ta)
            p = Any[]
            T0 = field(ds0, "t2m_ann", k0); T1 = field(ds, "t2m_ann", k)
            push!(p, ("ΔT annual [K]", T1 .- T0, (-25, 25), Reverse(:RdBu), identity))
            S0 = field(ds0, "t2m_sum", k0)
            S0 !== nothing && push!(p, ("ΔT summer [K]", field(ds, "t2m_sum", k) .- S0, (-25, 25), Reverse(:RdBu), identity))
            P0 = field(ds0, "pr_ann", k0)
            P0 !== nothing && push!(p, ("Precipitation ratio", field(ds, "pr_ann", k) ./ P0, (0.2, 5), :BrBG, log10))
            M0 = field(ds0, "smb_ann", k0)
            M0 !== nothing && push!(p, ("ΔSMB [m/a w.e.]", field(ds, "smb_ann", k) .- M0, (-3, 3), :RdBu, identity))
            rlab = isempty(opt["ref_run"]) ? "$(tlabel(ta)) minus $(tlabel(opt["ref"]))" :
                   "minus $(basename(normpath(opt["ref_run"])))"
            push!(rows, (rlab, p))
        end
    end
    ds0 === ds || close(ds0)
    close(ds)

    ncol = maximum(length(r[2]) for r in rows)
    w = 300
    h = round(Int, w / aspect)
    fig = Figure(size=(ncol * w + 40, length(rows) * (h + 110) + 40), fontsize=13)
    isempty(opt["title"]) || Label(fig[0, 1:ncol], opt["title"], fontsize=15, font=:bold)
    for (r, (rlab, panels)) in enumerate(rows)
        for (c, (ttl, v, cr, cm, sc)) in enumerate(panels)
            ax = Axis(fig[2r-1, c]; title=isempty(rlab) ? ttl : "$rlab\n$ttl", aspect=DataAspect(),
                      xticksvisible=false, yticksvisible=false,
                      xticklabelsvisible=false, yticklabelsvisible=false)
            hm = heatmap!(ax, x, y, v; colorrange=cr, colormap=cm, colorscale=sc,
                          nan_color=:transparent)
            cst === nothing || contour!(ax, cst...; levels=[0.5], color=:black, linewidth=0.6)
            Colorbar(fig[2r, c], hm; vertical=false, flipaxis=false)
        end
    end
    save(opt["out"], fig)
    println(opt["out"])
end

main(parse_args(ARGS))

# Focused analysis of the Yelmo-32KM ensemble: what does refining the
# marine_shelf grid do, relative to the m32 (mshlf-on-Yelmo-grid) baseline?
#
# Reads ONLY the JLD2 bundle written by analysis/export_mg_jld2.jl, so it runs
# locally with just JLD2 + CairoMakie -- no NetCDF, no access to the raw output.
#
# Reference   : y32KM_m32KM  (marine_shelf on the Yelmo grid = no refinement)
# Refinements : y32KM_m16KM, y32KM_m8KM  (finer marine_shelf grid)
# All three share the ANT-32KM Yelmo grid, so differences are point-by-point.
#
# Produces:
#   (1) mg_y32_refine_maps.png -- rows = bmb_shlf and H-H_obs; col 1 the m32
#       reference field, cols 2-3 the (refined - m32) difference;
#   (2) mg_y32_refine_ts.png   -- V_sle(t), mean shelf bmb(t), and a bmb_shlf
#       scatter (refined vs m32) with bias/RMSE.
# Also prints a scalar summary (ΔV_sle, RMSE of the bmb_shlf / H differences).
#
# Usage:
#   julia --project=analysis analysis/mg_y32_refine.jl [bundle.jld2] [out_dir]
#     bundle.jld2  default analysis/mg_ensemble.jld2
#     out_dir      default analysis/figures

using JLD2, CairoMakie, Printf, Statistics

const REF     = "y32KM_m32KM"
const REFINED = ["y32KM_m16KM", "y32KM_m8KM"]
mlabel(name) = replace(name, "y32KM_" => "")     # "y32KM_m8KM" -> "m8KM"

# --- masks from a run's final-time geometry --------------------------------
getmap(run, v) = run["maps"][v]
function masks(run)
    H  = getmap(run, "H_ice"); zs = getmap(run, "z_srf"); zb = getmap(run, "z_bed")
    ice  = H .> 1
    flt  = ice .& ((zs .- H) .> (zb .+ 1))
    grnd = ice .& .!flt
    return ice, flt, grnd
end

nanmask(f, m) = (g = Float64.(copy(f)); g[.!m] .= NaN; g)
function rmse_bias(d, m)
    v = d[m .& isfinite.(d)]
    isempty(v) && return (NaN, NaN, 0)
    return (sqrt(mean(v .^ 2)), mean(v), length(v))
end

# --- figure 1: reference field + refinement differences --------------------
function fig_maps(ens, out)
    runs = ens["runs"]
    ref  = runs[REF]
    x = ref["xc"]; y = ref["yc"]
    _, fltR, _ = masks(ref)
    obs = get(ens["obs"], ref["ygrid"], nothing)

    cols = 1 + length(REFINED)
    fig = Figure(size = (360 * cols, 720))

    # -- row 1: bmb_shlf (floating shelf only) --
    bR = nanmask(getmap(ref, "bmb_shlf"), fltR)
    ax = Axis(fig[1, 1]; aspect = DataAspect(), title = "bmb_shlf  $(mlabel(REF)) ref",
              ylabel = "yc [km]", titlesize = 11)
    hm = heatmap!(ax, x, y, bR; colormap = Reverse(:dense), colorrange = (-20, 0))
    Colorbar(fig[1, cols + 1], hm, label = "bmb_shlf [m/yr]")
    for (k, nm) in enumerate(REFINED)
        r = runs[nm]; _, flt, _ = masks(r)
        d = nanmask(getmap(r, "bmb_shlf") .- getmap(ref, "bmb_shlf"), flt .| fltR)
        rm, bi, n = rmse_bias(d, isfinite.(d))
        axd = Axis(fig[1, k + 1]; aspect = DataAspect(), titlesize = 10,
                   title = "Δbmb_shlf  $(mlabel(nm))-$(mlabel(REF))\nRMSE=$(round(rm,digits=2)) bias=$(round(bi,digits=2))")
        hmd = heatmap!(axd, x, y, d; colormap = :balance, colorrange = (-10, 10))
        k == length(REFINED) && Colorbar(fig[1, cols + 2], hmd, label = "Δ [m/yr]")
    end

    # -- row 2: H - H_obs (over ice), and ΔH for refinements --
    if obs !== nothing
        Hobs = Float64.(obs["H_ice"])
        iceR, _, _ = masks(ref)
        misR = nanmask(getmap(ref, "H_ice") .- Hobs, iceR .| (Hobs .> 10))
        rmH, biH, _ = rmse_bias(misR, isfinite.(misR))
        ax = Axis(fig[2, 1]; aspect = DataAspect(), ylabel = "yc [km]", xlabel = "xc [km]",
                  title = "H - H_obs  $(mlabel(REF)) ref\nRMSE=$(round(rmH,digits=1)) bias=$(round(biH,digits=1))", titlesize = 10)
        hm = heatmap!(ax, x, y, misR; colormap = :balance, colorrange = (-800, 800))
        Colorbar(fig[2, cols + 1], hm, label = "H - H_obs [m]")
    end
    for (k, nm) in enumerate(REFINED)
        r = runs[nm]; ice, _, _ = masks(r); iceR, _, _ = masks(ref)
        d = nanmask(getmap(r, "H_ice") .- getmap(ref, "H_ice"), ice .| iceR)
        rm, bi, _ = rmse_bias(d, isfinite.(d))
        axd = Axis(fig[2, k + 1]; aspect = DataAspect(), xlabel = "xc [km]", titlesize = 10,
                   title = "ΔH  $(mlabel(nm))-$(mlabel(REF))\nRMSE=$(round(rm,digits=1)) bias=$(round(bi,digits=1))")
        hmd = heatmap!(axd, x, y, d; colormap = :balance, colorrange = (-100, 100))
        k == length(REFINED) && Colorbar(fig[2, cols + 2], hmd, label = "ΔH [m]")
    end

    save(joinpath(out, "mg_y32_refine_maps.png"), fig)
end

# --- figure 2: timeseries + bmb_shlf scatter -------------------------------
function fig_ts(ens, out)
    runs = ens["runs"]
    names = [REF; REFINED]
    fig = Figure(size = (1200, 420))

    axV = Axis(fig[1, 1]; xlabel = "time [kyr]", ylabel = "V_sle [m]", title = "Ice volume", titlesize = 11)
    axB = Axis(fig[1, 2]; xlabel = "time [kyr]", ylabel = "bmb_shlf [m/yr]", title = "Mean shelf bmb", titlesize = 11)
    for nm in names
        ts = runs[nm]["ts"]
        haskey(ts, "V_sle")    && lines!(axV, ts["time"] ./ 1e3, ts["V_sle"];    label = mlabel(nm))
        haskey(ts, "bmb_shlf") && lines!(axB, ts["time"] ./ 1e3, ts["bmb_shlf"]; label = mlabel(nm))
    end
    axislegend(axV; position = :rb, labelsize = 9)

    # scatter: finest refinement bmb_shlf vs m32 (floating cells)
    ref = runs[REF]; _, fltR, _ = masks(ref)
    r8  = runs[REFINED[end]]; _, flt8, _ = masks(r8)
    m = fltR .& flt8
    xs = Float64.(getmap(ref, "bmb_shlf"))[m]; ys = Float64.(getmap(r8, "bmb_shlf"))[m]
    ax = Axis(fig[1, 3]; aspect = DataAspect(), titlesize = 10,
              xlabel = "bmb_shlf $(mlabel(REF)) [m/yr]", ylabel = "bmb_shlf $(mlabel(REFINED[end])) [m/yr]",
              title = "$(mlabel(REFINED[end])) vs $(mlabel(REF))  n=$(length(xs))")
    scatter!(ax, xs, ys; markersize = 4, color = (:steelblue, 0.4))
    lo, hi = -30.0, 2.0
    lines!(ax, [lo, hi], [lo, hi]; color = :black, linestyle = :dash)
    xlims!(ax, lo, hi); ylims!(ax, lo, hi)

    save(joinpath(out, "mg_y32_refine_ts.png"), fig)
end

# --- scalar summary --------------------------------------------------------
function print_summary(ens)
    runs = ens["runs"]; ref = runs[REF]
    vref = haskey(ref["ts"], "V_sle") ? ref["ts"]["V_sle"][end] : NaN
    @printf("%-12s  V_sle=%.3f m (reference)\n", mlabel(REF), vref)
    for nm in REFINED
        r = runs[nm]; _, flt, _ = masks(r); _, fltR, _ = masks(ref)
        db = (getmap(r, "bmb_shlf") .- getmap(ref, "bmb_shlf"))
        rb, bb, nb = rmse_bias(db, flt .| fltR)
        ice, _, _ = masks(r); iceR, _, _ = masks(ref)
        dh = (getmap(r, "H_ice") .- getmap(ref, "H_ice"))
        rh, bh, _ = rmse_bias(dh, ice .| iceR)
        v = haskey(r["ts"], "V_sle") ? r["ts"]["V_sle"][end] : NaN
        @printf("%-12s  ΔV_sle=%+.3f m   Δbmb RMSE=%.2f bias=%+.2f (n=%d)   ΔH RMSE=%.1f bias=%+.1f m\n",
                mlabel(nm), v - vref, rb, bb, nb, rh, bh)
    end
end

function main()
    bundle = length(ARGS) >= 1 ? ARGS[1] : "analysis/mg_ensemble.jld2"
    out    = length(ARGS) >= 2 ? ARGS[2] : "analysis/figures"
    mkpath(out)
    ens = load(bundle, "ens")
    fig_maps(ens, out)
    fig_ts(ens, out)
    print_summary(ens)
    println("wrote figures to $out")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

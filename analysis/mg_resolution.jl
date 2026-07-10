# Analysis for the Antarctica multigrid resolution test (yelmox).
#
# Reads the six runs staged by yelmox/run_mg_resolution.sh and produces:
#   (1) map of the shelf basal mass balance (yelmo%bnd%bmb_shlf) for the
#       reference run;
#   (2) comparison -- one row per completed run, left column the actual
#       bmb_shlf field, right column its difference vs. the reference
#       (bilinearly remapped onto the reference grid);
#   (3) permutation timeseries -- ice volume and mean shelf bmb over time, one
#       column per Yelmo grid, one line per marine_shelf resolution.
#
# Only runs that reached the 25 kyr end are used for the map/difference/scatter
# figures (run_complete); incomplete or crashed runs are skipped so their early
# slices are not mistaken for a spun-up state. The reference is the completed
# run with the FINEST marine_shelf grid (refrun) rather than a hard-wired
# ANT-8KM run, which may not have finished.
#
# All bmb_shlf comparisons use yelmo%bnd%bmb_shlf, i.e. the field aggregated
# onto the Yelmo grid (whatever the marine_shelf grid was).
#
# Usage (from the yelmox root):
#   julia --project=analysis analysis/mg_resolution.jl [run_root] [out_dir]
#     run_root  default output/mg_opt
#     out_dir   default analysis/figures

using CairoMakie, NCDatasets, Printf

# --- run matrix ------------------------------------------------------------
# name : output subdir; ygrid : Yelmo grid (=> filename suffix); mgrid : mshlf grid.
struct Run; name::String; ygrid::String; mgrid::String; end
runs = [
    Run("y8KM_m8KM",   "ANT-8KM",  "ANT-8KM"),
    Run("y16KM_m16KM", "ANT-16KM", "ANT-16KM"),
    Run("y32KM_m32KM", "ANT-32KM", "ANT-32KM"),
    Run("y32KM_m8KM",  "ANT-32KM", "ANT-8KM"),
    Run("y32KM_m16KM", "ANT-32KM", "ANT-16KM"),
    Run("y16KM_m8KM",  "ANT-16KM", "ANT-8KM"),
]
findrun(nm) = runs[findfirst(r -> r.name == nm, runs)]

# res label like "8KM" from a grid name / run
reslabel(grid) = replace(grid, "ANT-" => "")

# integer km of a grid name, e.g. "ANT-8KM" -> 8
gridkm(grid) = parse(Int, match(r"(\d+)KM", grid).captures[1])

# A run "completed" if its 2D file exists and reached (near) the 25 kyr end.
# Incomplete/crashed runs only wrote early slices; comparing their last slice
# against a spun-up run would mix model times, so they are skipped.
function run_complete(root, r; tmin = 24999.0)
    p = yelmo2d(root, r)
    isfile(p) || return false
    NCDataset(p) do ds
        haskey(ds, "time") || return false
        t = ds["time"][:]
        return length(t) > 0 && maximum(t) >= tmin
    end
end

# Reference run for the bmb_shlf comparison: the completed run with the FINEST
# marine_shelf grid. Was hard-wired to the ANT-8KM run, but that (like the other
# fine-Yelmo runs) may not have finished; this degrades gracefully to whatever
# finest-shelf run did complete. Returns nothing if no run completed.
function refrun(root)
    cand = filter(r -> run_complete(root, r), runs)
    isempty(cand) && return nothing
    return cand[argmin([gridkm(r.mgrid) for r in cand])]
end

# --- IO helpers ------------------------------------------------------------
yelmo2d(root, r) = joinpath(root, r.name, "yelmo.nc")
yelmots(root, r) = joinpath(root, r.name, "yelmo_ts.nc")

"Read (xc, yc, field[last time]) from a 2D file, or nothing if absent."
function read_map(path, var)
    isfile(path) || return nothing
    NCDataset(path) do ds
        haskey(ds, var) || return nothing
        x = Float64.(ds["xc"][:])
        y = Float64.(ds["yc"][:])
        f = Float64.(ds[var][:, :, end])
        # Geometry masks. bmb_shlf is only physically applied under floating ice;
        # over grounded ice / open ocean it is a non-physical phantom value. Also
        # build 0/1 masks of ice cover and grounded ice for PD outline contours.
        flt  = trues(size(f))
        ice  = zeros(size(f))       # ice-covered (grounded + floating)
        grnd = zeros(size(f))       # grounded ice
        if all(k -> haskey(ds, k), ("H_ice", "z_srf", "z_bed"))
            H  = Float64.(ds["H_ice"][:, :, end])
            zs = Float64.(ds["z_srf"][:, :, end])
            zb = Float64.(ds["z_bed"][:, :, end])
            icev = H .> 1
            flt  = icev .& ((zs .- H) .> (zb .+ 1))
            ice  = Float64.(icev)
            grnd = Float64.(icev .& .!flt)
        end
        return (x, y, f, flt, ice, grnd)
    end
end

# Overlay PD outlines: ice-cover margin + grounding line, as thin black contours.
function pd_outlines!(ax, x, y, ice, grnd)
    contour!(ax, x, y, ice;  levels = [0.5], color = :black, linewidth = 0.4)
    contour!(ax, x, y, grnd; levels = [0.5], color = :black, linewidth = 0.4)
end

"Read (time, var) from a timeseries file, or nothing if absent."
function read_ts(path, var)
    isfile(path) || return nothing
    NCDataset(path) do ds
        haskey(ds, var) || return nothing
        return (Float64.(ds["time"][:]), Float64.(ds[var][:]))
    end
end

# --- bilinear remap onto a target regular grid -----------------------------
"Bilinearly interpolate f (on ascending regular grid xs,ys) onto (xt,yt)."
function remap_bilin(xs, ys, f, xt, yt)
    out = fill(NaN, length(xt), length(yt))
    for (j, y) in enumerate(yt), (i, x) in enumerate(xt)
        (x < xs[1] || x > xs[end] || y < ys[1] || y > ys[end]) && continue
        ix = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
        iy = clamp(searchsortedlast(ys, y), 1, length(ys) - 1)
        tx = (x - xs[ix]) / (xs[ix+1] - xs[ix])
        ty = (y - ys[iy]) / (ys[iy+1] - ys[iy])
        out[i, j] = (1 - tx) * (1 - ty) * f[ix, iy]   + tx * (1 - ty) * f[ix+1, iy] +
                    (1 - tx) * ty       * f[ix, iy+1] + tx * ty       * f[ix+1, iy+1]
    end
    return out
end

# --- figure 1: reference bmb_shlf map --------------------------------------
function fig_ref_map(root, out)
    r = refrun(root)
    r === nothing && (@warn "no completed run, skipping fig 1"; return)
    m = read_map(yelmo2d(root, r), "bmb_shlf")
    m === nothing && (@warn "reference run missing, skipping fig 1"; return)
    x, y, f, flt, ice, grnd = m
    fm = copy(f); fm[.!flt] .= NaN                 # show floating shelf melt only

    fig = Figure(size = (620, 640))
    ax = Axis(fig[1, 1]; aspect = DataAspect(), title = "bmb_shlf  ($(r.name) reference)",
              xlabel = "xc [km]", ylabel = "yc [km]")
    hm = heatmap!(ax, x, y, fm; colormap = Reverse(:dense), colorrange = (-20, 0))
    pd_outlines!(ax, x, y, ice, grnd)
    Colorbar(fig[1, 2], hm, label = "bmb_shlf [m/yr]")
    save(joinpath(out, "mg_bmb_shlf_ref.png"), fig)
end

# --- figure 2: core comparison, actual + difference vs 8KM reference -------
# zoom = ((xmin,xmax),(ymin,ymax)) in km restricts the view (e.g. WAIS); tag is
# appended to the output filename.
function fig_core_compare(root, out; zoom = nothing, tag = "")
    ref = refrun(root)
    ref === nothing && (@warn "no completed run, skipping fig 2"; return)
    rm = read_map(yelmo2d(root, ref), "bmb_shlf")
    rm === nothing && (@warn "reference run missing, skipping fig 2"; return)
    xr, yr, fr, fltr, icer, grndr = rm
    fr = copy(fr); fr[.!fltr] .= NaN               # reference: floating shelf only

    # All completed runs (skip incomplete/crashed ones), compared vs reference.
    cores = filter(r -> run_complete(root, r), runs)

    setzoom!(ax) = zoom !== nothing && (xlims!(ax, zoom[1]...); ylims!(ax, zoom[2]...))

    fig = Figure(size = (860, 380 * length(cores)))
    for (row, r) in enumerate(cores)
        m = read_map(yelmo2d(root, r), "bmb_shlf")
        m === nothing && continue
        x, y, f, flt, ice, grnd = m
        f = copy(f); f[.!flt] .= NaN               # floating shelf melt only

        # left: actual field
        axl = Axis(fig[row, 1]; aspect = DataAspect(),
                   title = "bmb_shlf  ($(r.name))",
                   ylabel = "yc [km]")
        hml = heatmap!(axl, x, y, f; colormap = Reverse(:dense), colorrange = (-20, 0))
        pd_outlines!(axl, x, y, ice, grnd)
        setzoom!(axl)
        row == length(cores) && (axl.xlabel = "xc [km]")
        row == 1 && Colorbar(fig[1, 2], hml, label = "bmb_shlf [m/yr]")

        # right: difference vs reference (remap onto 8KM ref grid)
        fq = remap_bilin(x, y, f, xr, yr)
        d = fq .- fr
        axr = Axis(fig[row, 3]; aspect = DataAspect(),
                   title = "difference vs $(ref.name)")
        hmr = heatmap!(axr, xr, yr, d; colormap = :balance, colorrange = (-10, 10))
        pd_outlines!(axr, xr, yr, icer, grndr)     # reference geometry (diff is on ref grid)
        setzoom!(axr)
        row == length(cores) && (axr.xlabel = "xc [km]")
        row == 1 && Colorbar(fig[1, 4], hmr, label = "Δ bmb_shlf [m/yr]")
    end
    save(joinpath(out, "mg_bmb_shlf_core_compare$(tag).png"), fig)
end

# --- figure 3: permutation timeseries --------------------------------------
# One column per Yelmo grid; each line a marine_shelf resolution. Top row ice
# volume, bottom row mean shelf bmb.
function fig_perm_timeseries(root, out)
    groups = [("ANT-32KM", ["y32KM_m32KM", "y32KM_m16KM", "y32KM_m8KM"]),
              ("ANT-16KM", ["y16KM_m16KM", "y16KM_m8KM"])]

    fig = Figure(size = (1000, 700))
    for (col, (ygrid, names)) in enumerate(groups)
        axV = Axis(fig[1, col]; title = "Yelmo $(reslabel(ygrid))",
                   ylabel = col == 1 ? "V_ice [1e6 km³]" : "",
                   xlabel = "")
        axB = Axis(fig[2, col];
                   ylabel = col == 1 ? "bmb_shlf [m/yr]" : "",
                   xlabel = "time [yr]")
        nplot = 0
        for nm in names
            r = findrun(nm)
            lbl = "mshlf $(reslabel(r.mgrid))"
            dV = read_ts(yelmots(root, r), "V_ice")
            dV !== nothing && (lines!(axV, dV[1], dV[2]; label = lbl); nplot += 1)
            dB = read_ts(yelmots(root, r), "bmb_shlf")
            dB !== nothing && lines!(axB, dB[1], dB[2]; label = lbl)
        end
        nplot > 0 && axislegend(axV; position = :rt, labelsize = 10)
    end
    save(joinpath(out, "mg_perm_timeseries.png"), fig)
end

# --- area-mean aggregation of a fine field onto a coarser grid -------------
"Average fine field ff (grid xf,yf) onto coarse centers (xc,yc); NaN if empty."
function aggregate_mean(xf, yf, ff, xc, yc)
    ix = [argmin(abs.(xc .- v)) for v in xf]     # nearest coarse col per fine col
    iy = [argmin(abs.(yc .- v)) for v in yf]
    sums = zeros(length(xc), length(yc)); cnts = zeros(Int, length(xc), length(yc))
    for j in 1:length(yf), i in 1:length(xf)
        v = ff[i, j]; isnan(v) && continue
        sums[ix[i], iy[j]] += v; cnts[ix[i], iy[j]] += 1
    end
    out = fill(NaN, length(xc), length(yc))
    for k in eachindex(out); cnts[k] > 0 && (out[k] = sums[k] / cnts[k]); end
    return out
end

# --- figure 4: scatter of coarse vs 8KM reference (floating shelf) ----------
function fig_scatter(root, out)
    ref = refrun(root)
    ref === nothing && (@warn "no completed run, skipping scatter"; return)
    rm = read_map(yelmo2d(root, ref), "bmb_shlf")
    rm === nothing && (@warn "reference run missing, skipping scatter"; return)
    xr, yr, fr, fltr, _, _ = rm
    refm = copy(fr); refm[.!fltr] .= NaN           # reference floating shelf field

    others = filter(r -> run_complete(root, r) && r.name != ref.name, runs)
    isempty(others) && (@warn "no other completed run, skipping scatter"; return)
    lo, hi = -30.0, 2.0
    fig = Figure(size = (460 * length(others), 480))
    for (col, r) in enumerate(others)
        m = read_map(yelmo2d(root, r), "bmb_shlf")
        m === nothing && continue
        x, y, f, flt, _, _ = m
        fc = copy(f); fc[.!flt] .= NaN                       # coarse floating field
        ragg = aggregate_mean(xr, yr, refm, x, y)            # 8KM ref -> coarse (area mean)
        pair = .!isnan.(fc) .& .!isnan.(ragg)
        xs = ragg[pair]; ys = fc[pair]
        bias = sum(ys .- xs) / length(xs)
        rmse = sqrt(sum((ys .- xs) .^ 2) / length(xs))
        ax = Axis(fig[1, col]; aspect = DataAspect(),
                  title = "$(r.name) vs $(ref.name)   n=$(length(xs))  bias=$(round(bias, digits=2))  RMSE=$(round(rmse, digits=2))",
                  xlabel = "bmb_shlf  $(ref.name) ref (aggregated) [m/yr]",
                  ylabel = "bmb_shlf  $(r.name) [m/yr]", titlesize = 10)
        scatter!(ax, xs, ys; markersize = 4, color = (:steelblue, 0.4))
        lines!(ax, [lo, hi], [lo, hi]; color = :black, linestyle = :dash)
        xlims!(ax, lo, hi); ylims!(ax, lo, hi)
    end
    save(joinpath(out, "mg_bmb_shlf_scatter.png"), fig)
end

# --- main ------------------------------------------------------------------
function main()
    root = length(ARGS) >= 1 ? ARGS[1] : "output/mg_opt"
    out  = length(ARGS) >= 2 ? ARGS[2] : "analysis/figures"
    mkpath(out)
    fig_ref_map(root, out)
    fig_core_compare(root, out)
    fig_core_compare(root, out; zoom = ((-2308.0, 452.0), (-1732.0, 1476.0)), tag = "_WAIS")
    fig_scatter(root, out)
    fig_perm_timeseries(root, out)
    println("wrote figures to $out")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

# Analysis for the Antarctica multigrid resolution test (yelmox_mg).
#
# Reads the six runs staged by yelmox_mg/run_mg_resolution.sh and produces:
#   (1) map of the shelf basal mass balance (yelmo%bnd%bmb_shlf) for the
#       ANT-8KM-everywhere reference run;
#   (2) core comparison -- rows = 8/16/32 km core runs, left column the actual
#       bmb_shlf field, right column its difference vs. the ANT-8KM reference
#       (coarser runs bilinearly remapped onto the 8KM reference grid);
#   (3) permutation timeseries -- ice volume and mean shelf bmb over time, one
#       column per Yelmo grid, one line per marine_shelf resolution.
#
# All bmb_shlf comparisons use yelmo%bnd%bmb_shlf, i.e. the field aggregated
# onto the Yelmo grid (whatever the marine_shelf grid was).
#
# Usage (from the yelmox root):
#   julia --project=analysis analysis/mg_resolution.jl [run_root] [out_dir]
#     run_root  default output/mg
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

# --- IO helpers ------------------------------------------------------------
yelmo2d(root, r) = joinpath(root, r.name, "yelmo_$(r.ygrid).nc")
yelmots(root, r) = joinpath(root, r.name, "yelmo_$(r.ygrid)_ts.nc")

"Read (xc, yc, field[last time]) from a 2D file, or nothing if absent."
function read_map(path, var)
    isfile(path) || return nothing
    NCDataset(path) do ds
        haskey(ds, var) || return nothing
        x = Float64.(ds["xc"][:])
        y = Float64.(ds["yc"][:])
        f = Float64.(ds[var][:, :, end])
        return (x, y, f)
    end
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
    r = findrun("y8KM_m8KM")
    m = read_map(yelmo2d(root, r), "bmb_shlf")
    m === nothing && (@warn "reference run missing, skipping fig 1"; return)
    x, y, f = m
    fm = copy(f); fm[fm .== 0.0] .= NaN            # hide grounded / no-melt points

    fig = Figure(size = (620, 640))
    ax = Axis(fig[1, 1]; aspect = DataAspect(), title = "bmb_shlf  (ANT-8KM reference)",
              xlabel = "xc [km]", ylabel = "yc [km]")
    hm = heatmap!(ax, x, y, fm; colormap = :dense, colorrange = (-20, 0))
    Colorbar(fig[1, 2], hm, label = "bmb_shlf [m/yr]")
    save(joinpath(out, "mg_bmb_shlf_ref.png"), fig)
end

# --- figure 2: core comparison, actual + difference vs 8KM reference -------
function fig_core_compare(root, out)
    ref = findrun("y8KM_m8KM")
    rm = read_map(yelmo2d(root, ref), "bmb_shlf")
    rm === nothing && (@warn "reference run missing, skipping fig 2"; return)
    xr, yr, fr = rm

    cores = [findrun("y8KM_m8KM"), findrun("y16KM_m16KM"), findrun("y32KM_m32KM")]

    fig = Figure(size = (860, 1120))
    for (row, r) in enumerate(cores)
        m = read_map(yelmo2d(root, r), "bmb_shlf")
        m === nothing && continue
        x, y, f = m

        # left: actual field
        axl = Axis(fig[row, 1]; aspect = DataAspect(),
                   title = "bmb_shlf  ($(reslabel(r.ygrid)) core)",
                   ylabel = "yc [km]")
        fm = copy(f); fm[fm .== 0.0] .= NaN
        hml = heatmap!(axl, x, y, fm; colormap = :dense, colorrange = (-20, 0))
        row == length(cores) && (axl.xlabel = "xc [km]")
        row == 1 && Colorbar(fig[1, 2], hml, label = "bmb_shlf [m/yr]")

        # right: difference vs reference (remap onto 8KM ref grid)
        fq = remap_bilin(x, y, f, xr, yr)
        d = fq .- fr
        axr = Axis(fig[row, 3]; aspect = DataAspect(),
                   title = "difference vs ANT-8KM")
        hmr = heatmap!(axr, xr, yr, d; colormap = :balance, colorrange = (-10, 10))
        row == length(cores) && (axr.xlabel = "xc [km]")
        row == 1 && Colorbar(fig[1, 4], hmr, label = "Δ bmb_shlf [m/yr]")
    end
    save(joinpath(out, "mg_bmb_shlf_core_compare.png"), fig)
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

# --- main ------------------------------------------------------------------
function main()
    root = length(ARGS) >= 1 ? ARGS[1] : "output/mg"
    out  = length(ARGS) >= 2 ? ARGS[2] : "analysis/figures"
    mkpath(out)
    fig_ref_map(root, out)
    fig_core_compare(root, out)
    fig_perm_timeseries(root, out)
    println("wrote figures to $out")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

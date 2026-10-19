# Bias correction MAR - ERA5 on GreenlandPaleo, for [bias] of the snapshot
# backend (docs/_grl-paleo-plan.md, "Reference climate on GRL-PAL"): with ERA5
# as the present-day reference, the corrected climate is MAR over Greenland.
#
# Both datasets 1981-2010, monthly. The difference is taken at sea level, with
# the reduction of the snapshot backend (libs/snapshots.f90):
#   tsl = tas + zs*L(m),  pr_sl = pr*exp(f_p*zs*L(m)),
#   L(m) = l1 + (l1-l2)*cos(2*pi*(m-1)/12)   (north: l2 in July)
# so that
#   dtas(m) = tsl_MAR - tsl_ERA5            [K]
#   fpr(m)  = pr_sl_MAR / pr_sl_ERA5        [1], clamped to [fpr_min, fpr_max]
# are independent of elevation. Outside the MAR coverage (f_valid of t2m and pr
# >= 0.5) the bias of the nearest covered cell, smoothed over the coverage
# (gaussian of width sigma, log(fpr)), is used, weighted by w(d), a cosine from
# 1 at the coverage edge to 0 at the distance D:
#   dtas = w*dtas, fpr = 1 + w*(fpr - 1).
# The GRL grid is a crop of the GRL-PAL grid (same projection and resolution).
#
# Usage (from the yelmox root, fesmdata and ice_data links in place):
#   julia --project=scripts/grl-paleo scripts/grl-paleo/make_bias_mar_era5.jl \
#         [--res 16] [--l1 0.0080] [--l2 0.0065] [--f_p 0.06] [--D 100] [--sigma 50] \
#         [--fpr_min 0.2] [--fpr_max 5] [--fig]
#
# Writes ice_data/GreenlandPaleo/GRL-PAL-<res>KM/GRL-PAL-<res>KM_BIAS-MARv3.14-ERA5-1981-2010.nc
# (and a figure next to it with --fig). l1, l2, f_p must be the [snap] lapse and
# f_p of the par file; they are written as attributes.

using NCDatasets
using CairoMakie
using Printf

function parse_args(args)
    opt = Dict{String,Any}("res" => 16, "l1" => 0.0080, "l2" => 0.0065, "f_p" => 0.06,
                           "D" => 100.0, "sigma" => 50.0, "fpr_min" => 0.2, "fpr_max" => 5.0, "fig" => false)
    i = 1
    while i <= length(args)
        key = lstrip(args[i], '-')
        haskey(opt, key) || error("unknown option $(args[i])")
        if key == "fig"
            opt[key] = true
            i += 1
        else
            opt[key] = key == "res" ? parse(Int, args[i+1]) : parse(Float64, args[i+1])
            i += 2
        end
    end
    return opt
end

lapse_monthly(l1, l2) = [l1 + (l1 - l2) * cos(2pi * (m - 1) / 12) for m in 1:12]

# Read a variable as Float64 with missing -> NaN.
readvar(ds, name) = Float64.(coalesce.(Array(ds[name]), NaN))

function main(opt)
    res  = opt["res"]
    grl  = "GRL-$(res)KM"
    pal  = "GRL-PAL-$(res)KM"
    file_mar = "fesmdata/Greenland/$grl/$(grl)_SMB-MARv3.14-1981-2010.nc"
    file_era = "fesmdata/GreenlandPaleo/$pal/$(pal)_Atmosphere-ERA5-1981-2010.nc"
    dir_out  = "ice_data/GreenlandPaleo/$pal"
    file_out = joinpath(dir_out, "$(pal)_BIAS-MARv3.14-ERA5-1981-2010.nc")

    L   = lapse_monthly(opt["l1"], opt["l2"])
    f_p = opt["f_p"]

    # ERA5 on GRL-PAL: t2m [K], pr [kg m-2 d-1], zs [m]; arrays (xc, yc, month).
    dse = NCDataset(file_era)
    xc  = Float64.(dse["xc"][:]);  yc = Float64.(dse["yc"][:])
    te  = readvar(dse, "t2m");  pe = readvar(dse, "pr");  ze = readvar(dse, "zs")
    crs_attrib = Dict(dse["crs"].attrib)
    close(dse)

    # MAR on GRL: t2m [K], pr [kg m-2 yr-1, 1 yr = 365.25 d], z_srf [m].
    dsm = NCDataset(file_mar)
    xm  = Float64.(dsm["xc"][:]);  ym = Float64.(dsm["yc"][:])
    tm  = readvar(dsm, "t2m");  pm = readvar(dsm, "pr") ./ 365.25;  zm = readvar(dsm, "z_srf")
    fv  = min.(readvar(dsm, "f_valid_t2m"), readvar(dsm, "f_valid_pr"))
    close(dsm)

    # Index of the GRL crop in GRL-PAL.
    ix = [findfirst(x -> isapprox(x, v; atol=1e-3), xc) for v in xm]
    iy = [findfirst(y -> isapprox(y, v; atol=1e-3), yc) for v in ym]
    (any(isnothing, ix) || any(isnothing, iy)) && error("$grl is not a crop of $pal")
    ix = Int.(ix);  iy = Int.(iy)

    nx, ny = length(xc), length(yc)
    dtas0 = fill(NaN, nx, ny, 12)
    fpr0  = fill(NaN, nx, ny, 12)
    cov   = falses(nx, ny)

    for (jm, j) in enumerate(iy), (im, i) in enumerate(ix)
        fv[im, jm] >= 0.5 || continue
        cov[i, j] = true
        for m in 1:12
            tsl_m = tm[im, jm, m] + zm[im, jm] * L[m]
            tsl_e = te[i, j, m]   + ze[i, j]   * L[m]
            psl_m = pm[im, jm, m] * exp(f_p * zm[im, jm] * L[m])
            psl_e = pe[i, j, m]   * exp(f_p * ze[i, j]   * L[m])
            dtas0[i, j, m] = tsl_m - tsl_e
            fpr0[i, j, m]  = psl_e > 0 ? clamp(psl_m / psl_e, opt["fpr_min"], opt["fpr_max"]) :
                             (psl_m > 0 ? opt["fpr_max"] : 1.0)
        end
    end
    covm = repeat(cov, 1, 1, 12)
    (any(isnan, dtas0[covm]) || any(isnan, fpr0[covm])) && error("missing values in the covered cells")

    # The bias smoothed over the coverage, for the cells outside it.
    dtas_s = smooth_covered(dtas0, cov, xc, yc, opt["sigma"])
    fpr_s  = exp.(smooth_covered(log.(fpr0), cov, xc, yc, opt["sigma"]))

    # Nearest covered cell and its distance [km], the weight w(d).
    pts = [(i, j) for j in 1:ny, i in 1:nx if cov[i, j]]
    dist = zeros(nx, ny)
    w    = ones(nx, ny)
    dtas = similar(dtas0);  fpr = similar(fpr0)
    D = opt["D"]
    for j in 1:ny, i in 1:nx
        if cov[i, j]
            dtas[i, j, :] = dtas0[i, j, :];  fpr[i, j, :] = fpr0[i, j, :]
            continue
        end
        dmin = Inf;  kmin = (0, 0)
        for (ii, jj) in pts
            d = hypot(xc[ii] - xc[i], yc[jj] - yc[j])
            if d < dmin
                dmin = d;  kmin = (ii, jj)
            end
        end
        dist[i, j] = dmin
        w[i, j] = dmin >= D ? 0.0 : 0.5 * (1 + cos(pi * dmin / D))
        dtas[i, j, :] = w[i, j] .* dtas_s[kmin..., :]
        fpr[i, j, :]  = 1 .+ w[i, j] .* (fpr_s[kmin..., :] .- 1)
    end

    mkpath(dir_out)
    isfile(file_out) && rm(file_out)
    NCDataset(file_out, "c") do ds
        ds.attrib["title"]   = "Bias correction MAR v3.14 - ERA5, 1981-2010, on $pal"
        ds.attrib["source"]  = "yelmox scripts/grl-paleo/make_bias_mar_era5.jl"
        ds.attrib["input_mar"]  = basename(file_mar)
        ds.attrib["input_era5"] = basename(file_era)
        ds.attrib["method"]  = "difference at sea level: tsl = tas + zs*L(m), " *
                               "pr_sl = pr*exp(f_p*zs*L(m)), L(m) = l1 + (l1-l2)*cos(2*pi*(m-1)/12); " *
                               "outside the MAR coverage the nearest covered value times w(d), " *
                               "smoothed over the coverage (gaussian sigma, log fpr), " *
                               "w = cosine from 1 at the coverage edge to 0 at D"
        ds.attrib["lapse_l1"] = opt["l1"]
        ds.attrib["lapse_l2"] = opt["l2"]
        ds.attrib["f_p"]      = f_p
        ds.attrib["taper_D_km"] = D
        ds.attrib["smooth_sigma_km"] = opt["sigma"]
        ds.attrib["fpr_min"]  = opt["fpr_min"]
        ds.attrib["fpr_max"]  = opt["fpr_max"]
        ds.attrib["history"]  = "created $(Libc.strftime("%Y-%m-%d", time()))"

        defDim(ds, "xc", nx);  defDim(ds, "yc", ny);  defDim(ds, "month", 12)
        defVar(ds, "xc", xc, ("xc",), attrib=["units" => "km", "standard_name" => "projection_x_coordinate", "axis" => "X"])
        defVar(ds, "yc", yc, ("yc",), attrib=["units" => "km", "standard_name" => "projection_y_coordinate", "axis" => "Y"])
        defVar(ds, "month", Int32.(1:12), ("month",), attrib=["units" => "1", "long_name" => "month of the year"])
        defVar(ds, "crs", Int32(0), (), attrib=collect(crs_attrib))
        defVar(ds, "dtas", Float32.(dtas), ("xc", "yc", "month"),
               attrib=["units" => "K", "long_name" => "air temperature correction MAR - ERA5 (at sea level)",
                       "grid_mapping" => "crs"])
        defVar(ds, "fpr", Float32.(fpr), ("xc", "yc", "month"),
               attrib=["units" => "1", "long_name" => "precipitation factor MAR / ERA5 (at sea level)",
                       "grid_mapping" => "crs"])
        defVar(ds, "w", Float32.(w), ("xc", "yc"),
               attrib=["units" => "1", "long_name" => "weight of the bias (1 in the MAR coverage)",
                       "grid_mapping" => "crs"])
        defVar(ds, "dist", Float32.(dist), ("xc", "yc"),
               attrib=["units" => "km", "long_name" => "distance to the MAR coverage",
                       "grid_mapping" => "crs"])
    end
    @printf("%s\n  dtas [K]: %.2f .. %.2f (ann. mean over coverage %.2f)\n  fpr: %.2f .. %.2f\n",
            file_out, minimum(dtas), maximum(dtas), sum(dtas[covm]) / count(covm),
            minimum(fpr), maximum(fpr))

    opt["fig"] && plot_bias(file_out, xc, yc, dtas, fpr, w)
    return file_out
end

function plot_bias(file_out, xc, yc, dtas, fpr, w)
    fig = Figure(size=(1400, 600))
    panels = [("dtas January [K]", dtas[:, :, 1], (-10, 10), Reverse(:RdBu)),
              ("dtas July [K]",    dtas[:, :, 7], (-5, 5),   Reverse(:RdBu)),
              ("fpr annual mean",  vec_mean(fpr), (0.2, 5),  :BrBG),
              ("w",                w,             (0, 1),    :viridis)]
    for (k, (ttl, f, cr, cm)) in enumerate(panels)
        ax = Axis(fig[1, k]; title=ttl, aspect=DataAspect())
        hm = heatmap!(ax, xc, yc, f; colorrange=cr, colormap=cm,
                      colorscale=(k == 3 ? log10 : identity))
        Colorbar(fig[2, k], hm; vertical=false)
    end
    png = replace(file_out, ".nc" => ".png")
    save(png, fig)
    println("  figure: ", png)
end

# Gaussian smoothing (width sigma [km], cut at 3 sigma) of the covered cells of
# f (x, y, month), normalised by the covered weight; NaN outside the coverage.
function smooth_covered(f, cov, xc, yc, sigma)
    nx, ny, nm = size(f)
    dx = abs(xc[2] - xc[1])
    r  = ceil(Int, 3sigma / dx)
    k  = [exp(-((a*dx)^2 + (b*dx)^2) / (2sigma^2)) for a in -r:r, b in -r:r]
    out = fill(NaN, nx, ny, nm)
    for j in 1:ny, i in 1:nx
        cov[i, j] || continue
        wsum = 0.0
        acc  = zeros(nm)
        for b in -r:r, a in -r:r
            ii, jj = i + a, j + b
            (1 <= ii <= nx && 1 <= jj <= ny && cov[ii, jj]) || continue
            wk = k[a+r+1, b+r+1]
            wsum += wk
            acc .+= wk .* f[ii, jj, :]
        end
        out[i, j, :] = acc ./ wsum
    end
    return out
end

vec_mean(a) = dropdims(sum(a, dims=3), dims=3) ./ size(a, 3)

main(parse_args(ARGS))

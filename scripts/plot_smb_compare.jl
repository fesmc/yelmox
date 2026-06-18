#!/usr/bin/env julia
# Compare smb_simple accumulation-modifier variants (baseline/dist/lat/both).
# Expects four run directories under <base>/, each with a yelmo2D.nc, produced by
# runme with smb_simple and (optionally) the smb_simple.k_acc_* / facc_*_min
# overrides. Writes <base>/compare_smb.png.
#
# Usage:
#   julia scripts/plot_smb_compare.jl [base_dir]
# where base_dir defaults to <repo>/tmp/smb_accum_test.
#
# Uses the global env packages CairoMakie + NCDatasets.

using CairoMakie
using NCDatasets

const BASE = length(ARGS) >= 1 ? ARGS[1] :
             joinpath(dirname(@__DIR__), "tmp", "smb_accum_test")
const VARIANTS = ["baseline", "dist", "lat", "both"]
const TITLES = Dict(
    "baseline" => "baseline\n(flat c0)",
    "dist"     => "dist\nk_acc_dist=5e-4, floor 0.5",
    "lat"      => "lat\nk_acc_lat=0.02, phi_ref=50, floor 0.3",
    "both"     => "both",
)

# Load last time slice. NCDatasets returns smb as (xc, yc, time).
function load_run(v)
    ds = NCDataset(joinpath(BASE, v, "yelmo2D.nc"))
    smb = Float64.(ds["smb"][:, :, end])     # m ie/yr
    H   = Float64.(ds["H_ice"][:, :, end])
    xc  = Float64.(ds["xc"][:])
    yc  = Float64.(ds["yc"][:])
    close(ds)
    return smb, H, xc, yc
end

# Accumulation only (smb>0 -> mm ie/yr; else NaN)
acc(s) = map(x -> x > 0 ? x * 1e3 : NaN, s)

function main()
    smb = Dict(v => load_run(v)[1] for v in VARIANTS)
    _, H, xc, yc = load_run("baseline")

    accv = Dict(v => acc(smb[v]) for v in VARIANTS)
    vmax = maximum(maximum(filter(!isnan, accv[v])) for v in VARIANTS)

    b = accv["baseline"]
    diffs = Dict(v => accv[v] .- b for v in ["dist", "lat", "both"])
    dmax = maximum(maximum(abs.(filter(!isnan, diffs[v]))) for v in keys(diffs))

    fig = Figure(size = (1600, 880))
    Label(fig[0, 1:4],
          "smb_simple accumulation realism: LIS-32KM, 10 yr, ice sheet & isostasy off " *
          "(contour = reference ice extent)";
          fontsize = 16, font = :bold)

    noticks = (xticksvisible = false, yticksvisible = false,
               xticklabelsvisible = false, yticklabelsvisible = false)

    # Row 1: accumulation field per variant
    hm = nothing
    for (j, v) in enumerate(VARIANTS)
        ax = Axis(fig[1, j]; title = TITLES[v], aspect = DataAspect(), noticks...)
        hm = heatmap!(ax, xc, yc, accv[v]; colorrange = (0, vmax), colormap = :viridis)
        contour!(ax, xc, yc, H; levels = [1.0], color = :white, linewidth = 0.7)
    end
    Colorbar(fig[1, 5], hm; label = "accumulation (mm ie/yr)")

    # Row 2: difference vs baseline (col 1 = caption)
    Label(fig[2, 1], "difference\nvs baseline\n(accum, mm ie/yr)";
          fontsize = 14, tellheight = false)
    hd = nothing
    for (j, v) in enumerate(["dist", "lat", "both"])
        ax = Axis(fig[2, j + 1]; title = "$v − baseline", aspect = DataAspect(), noticks...)
        hd = heatmap!(ax, xc, yc, diffs[v]; colorrange = (-dmax, dmax), colormap = :RdBu)
        contour!(ax, xc, yc, H; levels = [1.0], color = :black, linewidth = 0.7)
    end
    Colorbar(fig[2, 5], hd; label = "Δ accumulation (mm ie/yr)")

    out = joinpath(BASE, "compare_smb.png")
    save(out, fig)
    println("wrote ", out)
end

main()

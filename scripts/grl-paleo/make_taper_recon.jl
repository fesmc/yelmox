# Taper weight of the Greenland climate reconstructions on GreenlandPaleo, for
# [taper] (docs/_grl-paleo-plan.md, decision 5): now = ref + w*(now - ref).
# The reconstructions (Buizert et al., 2018; Badgeley et al., 2020) are on the
# Greenland grid only; w = 1 where all of them cover a cell (f_valid >= 0.5),
# then a cosine from 1 at the coverage edge to 0 at the distance D.
# The GRL grid is a crop of the GRL-PAL grid (same projection and resolution).
#
# Usage (from the yelmox root, fesmdata and ice_data links in place):
#   julia --project=scripts/grl-paleo scripts/grl-paleo/make_taper_recon.jl [--res 16] [--D 100] [--fig]
#
# Writes ice_data/GreenlandPaleo/GRL-PAL-<res>KM/GRL-PAL-<res>KM_TAPER-Buizert2018-Badgeley2020.nc
# (and a figure next to it with --fig).

using NCDatasets
using CairoMakie
using Printf

const PRODUCTS = ["Paleoclimate-Buizert2018", "Paleoclimate-Badgeley2020"]

function parse_args(args)
    opt = Dict{String,Any}("res" => 16, "D" => 100.0, "fig" => false)
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

readvar(ds, name) = Float64.(coalesce.(Array(ds[name]), NaN))

function main(opt)
    res = opt["res"]
    grl = "GRL-$(res)KM"
    pal = "GRL-PAL-$(res)KM"
    file_grid = "fesmdata/GreenlandPaleo/$pal/$(pal)_grid.nc"
    dir_out   = "ice_data/GreenlandPaleo/$pal"
    file_out  = joinpath(dir_out, "$(pal)_TAPER-Buizert2018-Badgeley2020.nc")

    dsg = NCDataset(file_grid)
    xc = Float64.(dsg["xc"][:]);  yc = Float64.(dsg["yc"][:])
    crs_attrib = Dict(dsg["crs"].attrib)
    close(dsg)
    nx, ny = length(xc), length(yc)

    # Coverage of all products on GRL-PAL (false outside the GRL crop).
    cov = falses(nx, ny)
    first = true
    for p in PRODUCTS
        ds = NCDataset("fesmdata/Greenland/$grl/$(grl)_$p.nc")
        xm = Float64.(ds["xc"][:]);  ym = Float64.(ds["yc"][:])
        fv = readvar(ds, "f_valid")
        close(ds)
        ix = [findfirst(x -> isapprox(x, v; atol=1e-3), xc) for v in xm]
        iy = [findfirst(y -> isapprox(y, v; atol=1e-3), yc) for v in ym]
        (any(isnothing, ix) || any(isnothing, iy)) && error("$grl is not a crop of $pal")
        c = falses(nx, ny)
        c[Int.(ix), Int.(iy)] .= fv .>= 0.5
        cov = first ? c : cov .& c
        first = false
    end

    # Distance [km] to the coverage, the weight w(d).
    pts = [(xc[i], yc[j]) for j in 1:ny, i in 1:nx if cov[i, j]]
    D = opt["D"]
    dist = zeros(nx, ny)
    w    = ones(nx, ny)
    for j in 1:ny, i in 1:nx
        cov[i, j] && continue
        d = minimum(p -> hypot(p[1] - xc[i], p[2] - yc[j]), pts)
        dist[i, j] = d
        w[i, j] = d >= D ? 0.0 : 0.5 * (1 + cos(pi * d / D))
    end

    mkpath(dir_out)
    isfile(file_out) && rm(file_out)
    NCDataset(file_out, "c") do ds
        ds.attrib["title"]   = "Taper weight of the Greenland reconstructions (Buizert 2018, Badgeley 2020) on $pal"
        ds.attrib["source"]  = "yelmox scripts/grl-paleo/make_taper_recon.jl"
        ds.attrib["inputs"]  = join(["$(grl)_$p.nc" for p in PRODUCTS], ", ")
        ds.attrib["method"]  = "w = 1 where all products cover a cell (f_valid >= 0.5), " *
                               "cosine from 1 at the coverage edge to 0 at D"
        ds.attrib["taper_D_km"] = D
        ds.attrib["history"] = "created $(Libc.strftime("%Y-%m-%d", time()))"

        defDim(ds, "xc", nx);  defDim(ds, "yc", ny)
        defVar(ds, "xc", xc, ("xc",), attrib=["units" => "km", "standard_name" => "projection_x_coordinate", "axis" => "X"])
        defVar(ds, "yc", yc, ("yc",), attrib=["units" => "km", "standard_name" => "projection_y_coordinate", "axis" => "Y"])
        defVar(ds, "crs", Int32(0), (), attrib=collect(crs_attrib))
        defVar(ds, "w", Float32.(w), ("xc", "yc"),
               attrib=["units" => "1", "long_name" => "weight of the reconstructions (1 in their coverage)",
                       "grid_mapping" => "crs"])
        defVar(ds, "dist", Float32.(dist), ("xc", "yc"),
               attrib=["units" => "km", "long_name" => "distance to the coverage of the reconstructions",
                       "grid_mapping" => "crs"])
    end
    @printf("%s\n  covered cells: %d of %d, w > 0: %d\n", file_out, count(cov), nx*ny, count(w .> 0))

    if opt["fig"]
        fig = Figure(size=(500, 600))
        ax = Axis(fig[1, 1]; title="w", aspect=DataAspect())
        hm = heatmap!(ax, xc, yc, w; colorrange=(0, 1), colormap=:viridis)
        Colorbar(fig[1, 2], hm)
        png = replace(file_out, ".nc" => ".png")
        save(png, fig)
        println("  figure: ", png)
    end
    return file_out
end

main(parse_args(ARGS))

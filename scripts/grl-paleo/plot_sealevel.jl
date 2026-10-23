# Global sea level: the FesmData SeaLevel products, the yelmox input/*.dat curves
# (dashed) and the curve of the GRL paleo setups, input/sealevel_spratt2016.dat
# (black; make_sealevel_spratt2016.jl).
#
# Usage (from the yelmox root, fesmdata link in place):
#   julia --project=scripts/grl-paleo scripts/grl-paleo/plot_sealevel.jl docs/figures/grlpal_sealevel.png
using NCDatasets, CairoMakie, DelimitedFiles
const D = "fesmdata/Records/SeaLevel"
recs = [("Spratt2016", "sl_short"), ("Waelbroeck2002", "sl"), ("Bintanja2008", "sl"), ("deBoer2014", "sl"),
        ("Berends2021", "sl"), ("Rohling2021", "sl"), ("Lambeck2014", "sl"), ("Grant2014-probabilistic", "rsl"),
        ("PaleoMIST", "sl")]
function rec(p, v)
    ds = NCDataset(joinpath(D, "SeaLevel-$p.nc"))
    t = -Float64.(ds["time"][:]) ./ 1000
    s = Float64.(coalesce.(vec(ds[v][:, :]), NaN))
    close(ds)
    i = sortperm(t); return t[i], s[i]
end
function dat(f)
    d = readdlm(f; comments = true, comment_char = '#', skipstart = startswith(readline(f), "time") ? 1 : 0)
    return Float64.(d[:, 1]) ./ 1000, Float64.(d[:, 2])
end
fig = Figure(size = (1300, 1000), fontsize = 14)
axs = [Axis(fig[1, 1:2]; title = "Global sea level, 0-150 ka", ylabel = "sea level [m]", xlabel = "time [ka]"),
       Axis(fig[2, 1]; title = "0-35 ka", ylabel = "sea level [m]", xlabel = "time [ka]"),
       Axis(fig[2, 2]; title = "LGM, 15-30 ka", ylabel = "sea level [m]", xlabel = "time [ka]")]
cols = Makie.wong_colors()
pal = vcat(cols, [:teal, :gray50, :olive, :brown])
for (k, (p, v)) in enumerate(recs)
    t, s = rec(p, v)
    for ax in axs
        lines!(ax, t, s; color = pal[k], label = "FesmData $p", linewidth = 1.2)
    end
end
for (k, (f, lab)) in enumerate([("input/sealevel_rohling_450kyr.dat", "yelmox rohling_450kyr (A)"),
                                ("input/sl_deglaciation.dat", "yelmox sl_deglaciation (B)"),
                                ("input/sealevel_waelbroeck2002.dat", "yelmox waelbroeck2002")])
    t, s = dat(f)
    for ax in axs
        lines!(ax, t, s; color = [:red, :magenta, :darkorange][k], linewidth = 3, linestyle = :dash, label = lab)
    end
end
t, sl = dat("input/sealevel_spratt2016.dat")
for ax in axs
    lines!(ax, t, sl; color = :black, linewidth = 3.5, label = "GRL paleo: Spratt2016, Holocene offset tapered")
end
xlims!(axs[1], -150, 0); xlims!(axs[2], -35, 0); xlims!(axs[3], -30, -15)
ylims!(axs[1], -160, 20); ylims!(axs[2], -160, 20); ylims!(axs[3], -160, -80)
Legend(fig[1:2, 3], axs[1]; framevisible = false, labelsize = 12)
save(ARGS[1], fig)
println(ARGS[1])

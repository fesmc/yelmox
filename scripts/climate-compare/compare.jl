# Compare the climate 2D files of two climate-only runs (A and B, from run.sh):
# per time, the domain means of the annual air temperature (t2m_ann), the annual
# precipitation (pr_ann) and the annual ocean temperature (to_ann, depth mean or
# one level), over all cells and over the cells with ice; maps of A, B and B-A at
# the first, middle and last times, and the time series of the means.
#
# Usage:
#   julia --project=scripts/climate-compare scripts/climate-compare/compare.jl OUTDIR \
#         [--file snap.nc] [--depth D] [--ext png]
#
#   OUTDIR    folder with the runs A/ and B/; the figures are written here
#   --file    the climate file in A/ and B/ (default: the first of snap.nc,
#             esm.nc, rembo.nc found in A/)
#   --depth   ocean level [m] for to_ann (nearest level; default: depth mean)
#   --ext     figure format, png or pdf (default png)
#
# The ice mask is H_ice > 0 from the climate file when it has H_ice, else from
# the first record of A/yelmo.nc when it is on the same grid; without either,
# the "ice" means are over all cells.

using NCDatasets
using CairoMakie
using Printf
using Statistics

const FIELDS = [
    (name = "t2m_ann", label = "Air temperature (ann)", units = "K",    cmap = :thermal),
    (name = "pr_ann",  label = "Precipitation (ann)",   units = "mm/a", cmap = :YlGnBu),
    (name = "to_ann",  label = "Ocean temperature (ann)", units = "K",  cmap = :thermal),
]

# ---- command line ------------------------------------------------------------

function parse_args(args)
    opts = Dict{String,String}()
    pos = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if startswith(a, "--")
            opts[a[3:end]] = args[i+1]
            i += 2
        else
            push!(pos, a)
            i += 1
        end
    end
    length(pos) == 1 || error("usage: compare.jl OUTDIR [--file snap.nc] [--depth D] [--ext png]")
    return pos[1], opts
end

function find_clim_file(dir)
    for f in ("snap.nc", "esm.nc", "rembo.nc")
        isfile(joinpath(dir, f)) && return f
    end
    error("no climate file (snap.nc, esm.nc, rembo.nc) in $dir")
end

# ---- reading -----------------------------------------------------------------

# A variable as Float64 with missing values as NaN.
readvar(ds, name) = Float64.(coalesce.(Array(ds[name]), NaN))

# One run: the times and, per field present in the file, an array (nx, ny, nt);
# to_ann reduced over depth (thickness-weighted mean, or the level nearest `depth`).
function read_run(path, depth)
    ds = NCDataset(path)
    run = Dict{String,Any}()
    run["time"] = readvar(ds, "time")
    run["xc"]   = readvar(ds, "xc")
    run["yc"]   = readvar(ds, "yc")
    for f in FIELDS
        haskey(ds, f.name) || continue
        v = readvar(ds, f.name)
        if f.name == "to_ann"
            v = reduce_depth(v, readvar(ds, "depth"), depth)
        end
        run[f.name] = v
    end
    run["H_ice"] = haskey(ds, "H_ice") ? readvar(ds, "H_ice")[:, :, 1] : nothing
    close(ds)
    return run
end

function reduce_depth(v, z, depth)
    if depth === nothing
        # Layer thickness around each level (midpoints, ends at the first/last level).
        zm = [z[1]; (z[1:end-1] .+ z[2:end]) ./ 2; z[end]]
        w  = diff(zm)
        w  = w ./ sum(w)
        return dropdims(sum(v .* reshape(w, 1, 1, :, 1), dims = 3), dims = 3)
    else
        k = argmin(abs.(z .- depth))
        @printf("to_ann at the level %.1f m (asked %.1f m)\n", z[k], depth)
        return v[:, :, k, :]
    end
end

function ice_mask(runA, dirA)
    H = runA["H_ice"]
    if H === nothing
        fnm = joinpath(dirA, "yelmo.nc")
        if isfile(fnm)
            NCDataset(fnm) do ds
                if haskey(ds, "H_ice")
                    H = readvar(ds, "H_ice")[:, :, 1]
                end
            end
        end
    end
    nx, ny = length(runA["xc"]), length(runA["yc"])
    if H === nothing || size(H) != (nx, ny)
        println("No ice mask on the climate grid: the ice means are over all cells.")
        return trues(nx, ny), false
    end
    return H .> 0, true
end

# ---- statistics --------------------------------------------------------------

nanmean(x) = (y = filter(!isnan, vec(x)); isempty(y) ? NaN : mean(y))

# Domain means per time over the cells in `mask`.
domain_means(v, mask) = [nanmean(view(v, :, :, k)[mask]) for k in axes(v, 3)]

function print_table(f, t, mA, mB, mAi, mBi, has_mask)
    println()
    println("$(f.name) [$(f.units)]  $(f.label)")
    ice = has_mask ? "ice cells" : "all cells (no mask)"
    @printf("%12s | %10s %10s %10s | %10s %10s %10s\n", "", "all cells", "", "", ice, "", "")
    @printf("%12s | %10s %10s %10s | %10s %10s %10s\n", "time [yr]", "A", "B", "B-A", "A", "B", "B-A")
    println("-"^82)
    for k in eachindex(t)
        @printf("%12.0f | %10.3f %10.3f %10.3f | %10.3f %10.3f %10.3f\n",
                t[k], mA[k], mB[k], mB[k] - mA[k], mAi[k], mBi[k], mBi[k] - mAi[k])
    end
end

# ---- plotting ----------------------------------------------------------------

function symmetric_range(x)
    m = maximum(abs, filter(!isnan, vec(x)); init = 0.0)
    return m > 0 ? (-m, m) : (-1.0, 1.0)
end

function value_range(a, b)
    v = filter(!isnan, [vec(a); vec(b)])
    isempty(v) && return (0.0, 1.0)
    lo, hi = quantile(v, (0.01, 0.99))
    return lo < hi ? (lo, hi) : (lo - 1, hi + 1)
end

function plot_field(f, t, xc, yc, vA, vB, mA, mB, mAi, mBi, has_mask, fnm)
    nt = length(t)
    ks = unique([1, (nt + 1) ÷ 2, nt])
    crange = value_range(vA[:, :, ks], vB[:, :, ks])
    drange = symmetric_range(vB[:, :, ks] .- vA[:, :, ks])

    fig = Figure(size = (1000, 330 * length(ks) + 350))
    Label(fig[0, 1:3], "$(f.label) [$(f.units)]", fontsize = 18, font = :bold)

    for (r, k) in enumerate(ks)
        tl = @sprintf("%.0f yr", t[k])
        for (c, (v, ttl)) in enumerate(((vA[:, :, k], "A"), (vB[:, :, k], "B"),
                                        (vB[:, :, k] .- vA[:, :, k], "B - A")))
            ax = Axis(fig[r, c], title = "$ttl, $tl", aspect = DataAspect())
            hidedecorations!(ax)
            if c < 3
                hm = heatmap!(ax, xc, yc, v, colormap = f.cmap, colorrange = crange,
                              nan_color = :grey85)
            else
                hm = heatmap!(ax, xc, yc, v, colormap = :RdBu, colorrange = drange,
                              nan_color = :grey85)
            end
            if r == length(ks) && c != 2
                Colorbar(fig[r+1, c], hm, vertical = false, flipaxis = false, width = Relative(0.8))
            end
        end
    end

    # Time series of the domain means.
    rt = length(ks) + 2
    reg = has_mask ? "ice" : "all"
    ax1 = Axis(fig[rt, 1:2], xlabel = "Time [yr]", ylabel = "Mean [$(f.units)]",
               title = "Domain mean")
    lines!(ax1, t, mA, color = :black, label = "A, all")
    lines!(ax1, t, mB, color = :red,   label = "B, all")
    if has_mask
        lines!(ax1, t, mAi, color = :black, linestyle = :dash, label = "A, $reg")
        lines!(ax1, t, mBi, color = :red,   linestyle = :dash, label = "B, $reg")
    end
    Legend(fig[rt+1, 1:3], ax1, orientation = :horizontal, framevisible = false)

    ax2 = Axis(fig[rt, 3], xlabel = "Time [yr]", ylabel = "B - A [$(f.units)]",
               title = has_mask ? "B - A of the means (dashed: ice)" : "B - A of the means")
    hlines!(ax2, [0.0], color = :grey60)
    lines!(ax2, t, mB .- mA, color = :blue)
    has_mask && lines!(ax2, t, mBi .- mAi, color = :blue, linestyle = :dash)

    save(fnm, fig)
    println("Wrote $fnm")
end

# ---- main --------------------------------------------------------------------

function main(args)
    outdir, opts = parse_args(args)
    dirA, dirB = joinpath(outdir, "A"), joinpath(outdir, "B")
    file  = get(opts, "file", find_clim_file(dirA))
    depth = haskey(opts, "depth") ? parse(Float64, opts["depth"]) : nothing
    ext   = get(opts, "ext", "png")

    runA = read_run(joinpath(dirA, file), depth)
    runB = read_run(joinpath(dirB, file), depth)
    mask, has_mask = ice_mask(runA, dirA)
    println("A: $(joinpath(dirA, file))\nB: $(joinpath(dirB, file))")
    has_mask && println("Ice mask: $(count(mask)) of $(length(mask)) cells")

    # The times in both files.
    t  = intersect(runA["time"], runB["time"])
    kA = [findfirst(==(x), runA["time"]) for x in t]
    kB = [findfirst(==(x), runB["time"]) for x in t]
    isempty(t) && error("A and B have no time in common")

    for f in FIELDS
        if !(haskey(runA, f.name) && haskey(runB, f.name))
            missing_in = join([l for (l, r) in (("A", runA), ("B", runB)) if !haskey(r, f.name)], ", ")
            println("\n$(f.name): not in $missing_in, skipped")
            continue
        end
        vA = runA[f.name][:, :, kA]
        vB = runB[f.name][:, :, kB]
        mA, mB   = domain_means(vA, trues(size(mask))), domain_means(vB, trues(size(mask)))
        mAi, mBi = domain_means(vA, mask), domain_means(vB, mask)
        print_table(f, t, mA, mB, mAi, mBi, has_mask)
        plot_field(f, t, runA["xc"], runA["yc"], vA, vB, mA, mB, mAi, mBi, has_mask,
                   joinpath(outdir, "climcmp_$(f.name).$ext"))
    end
end

main(ARGS)

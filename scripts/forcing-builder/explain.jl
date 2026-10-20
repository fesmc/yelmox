# explain.jl — what a run will impose, on the command line: the configuration
# of a base par file plus runme -p overrides (or of an existing run folder)
# resolved against the forcing library and the schema defaults, printed with
# the source of every forcing value, and the forcing drawn on the run's
# timeline (the forcing clock, each index, the homogeneous anomalies, the sea
# level). The same resolution and formulas as the forcing builder page
# (docs/forcing-builder.md).
#
# Usage (from the repository root):
#   julia --project=scripts/forcing-builder scripts/forcing-builder/explain.jl PAR [-p KEY=VALUE ...] [--fig OUT.png]
#   julia --project=scripts/forcing-builder scripts/forcing-builder/explain.jl RUNDIR [--fig OUT.png]
#
#   PAR     a par file, e.g. yelmox/yelmox_Greenland.nml
#   RUNDIR  a run folder: its staged par file (yelmox*.nml) is read, no overrides
#   -p      overrides as for runme: g.k=v a scalar, g.k==a,b a vector
#   --fig   the figure (default: forcing.png next to a run folder's par file,
#           else forcing-builder.png in the working directory); "none" for no figure

include(joinpath(@__DIR__, "Namelist.jl"))
include(joinpath(@__DIR__, "TsGen.jl"))
include(joinpath(@__DIR__, "Resolve.jl"))

using .Namelist, .TsGen, .Resolve
using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function parse_args(args)
    isempty(args) && error("usage: explain.jl PAR|RUNDIR [-p KEY=VALUE ...] [--fig OUT.png]")
    target = args[1]
    p = String[]; fig = nothing
    i = 2
    while i <= length(args)
        a = args[i]
        if a == "-p"
            i += 1
            while i <= length(args) && !startswith(args[i], "-")
                push!(p, args[i]); i += 1
            end
        elseif a == "--fig"
            fig = args[i+1]; i += 2
        else
            error("unknown argument $a")
        end
    end
    return target, p, fig
end

val(groups, g, k, d=nothing) = haskey(groups, g) ? value_of(groups[g], k, d) : d
names_of(v) = v isa AbstractVector ? [string(x) for x in v if !isempty(string(x))] : (isempty(string(v)) ? String[] : [string(v)])
fnum(x) = x isa Integer ? string(x) : (x isa Real ? @sprintf("%g", x) : string(x))

"The forcing clock: the timeline group, the model times and the forcing times of the run."
function timeline(groups)
    T = string(val(groups, "ctrl", "run_step", "ctrl"))
    haskey(groups, T) || (T = "ctrl")
    method = string(val(groups, T, "tstep_method", "const"))
    t0 = Float64(val(groups, T, "time_init", 0.0)); t1 = Float64(val(groups, T, "time_end", 0.0))
    tc = Float64(val(groups, T, "tstep_const", 0.0))
    cref = Float64(val(groups, "ctrl", "calendar_ref", 1950.0))
    t = range(t0, t1, length=401)
    tf = method == "const" ? fill(tc, length(t)) : (method == "cal" ? collect(t) .- cref : collect(t))
    return (T=T, method=method, t0=t0, t1=t1, tc=tc, cref=cref, t=collect(t), tf=tf)
end

"The index group as a Dict for TsGen, with norm active when the group itself sets norm_t0 and norm_t1."
function index_par(rg::RGroup)
    p = Dict{String,Any}(k => e.value for (k, e) in rg.entries)
    p["norm"] = all(k -> haskey(rg.entries, k) && rg.entries[k].source != "default", ("norm_t0", "norm_t1"))
    return p
end

function load_index_series(p, dirs)
    p["method"] == "series" || return nothing
    file = string(get(p, "series_file", ""))
    isempty(file) && return nothing
    for d in dirs
        path = joinpath(d, file)
        isfile(path) || continue
        ser = load_series(path; var=string(get(p, "series_var", "")), time_name=string(get(p, "series_time_name", "time")))
        p["norm"] && (ser = normalize(ser, p["norm_t0"], p["norm_t1"]))
        return ser
    end
    @warn "series $file not found"
    return nothing
end

function describe_index(io, name, rg, p, ser, tl, f)
    m = string(p["method"])
    println(io, "  [idx_$name] ($(rg.source)) method = $m")
    if m == "series"
        println(io, "    f = S(t_f), S: ", p["series_file"], p["norm"] ? @sprintf("; normalised: 0 at %s, 1 at %s", string(p["norm_t0"]), string(p["norm_t1"])) : "")
    elseif is_feedback(m)
        @printf(io, "    response-driven (%s): f in [%s, %s], |df/dt| <= %s, tol = %s%s\n", m, fnum(p["f_min"]), fnum(p["f_max"]), fnum(p["df_dt_max"]), fnum(p["tol"]), get(p, "with_kill", false) ? ", with_kill" : "")
    else
        s = Float64(p["df_sign"]) >= 0 ? 1 : -1
        @printf(io, "    f_min = %s, f_max = %s, df_sign = %d, dt_init = %s", fnum(p["f_min"]), fnum(p["f_max"]), s, fnum(p["dt_init"]))
        m in ("const", "ramp-slope") && @printf(io, ", df_dt_max = %s /yr", fnum(p["df_dt_max"]))
        m in ("ramp-time", "ramp-time-step", "sin") && @printf(io, ", dt_ramp = %s yr", fnum(p["dt_ramp"]))
        m == "ramp-time-step" && @printf(io, ", dt_conv = %s yr, f_conv = %s", fnum(p["dt_conv"]), fnum(p["f_conv"]))
        println(io)
        tl.method == "const" && println(io, "    (the forcing clock is constant: the index stays at its start value)")
    end
    Float64(get(p, "sigma", 0.0)) > 0 && println(io, "    noise sigma = ", fnum(p["sigma"]))
    Float64(get(p, "win_t1", 0.0)) > Float64(get(p, "win_t0", 0.0)) && @printf(io, "    window %s < t_f < %s: f -> %s + %s (f - %s)\n", fnum(p["win_t0"]), fnum(p["win_t1"]), fnum(p["win_ref"]), fnum(p["win_fac"]), fnum(p["win_ref"]))
    fin = filter(isfinite, f)
    isempty(fin) || @printf(io, "    over the run: f = %.4g at the start, %.4g at the end, range [%.4g, %.4g]\n", f[1] + 0.0, f[end] + 0.0, minimum(fin), maximum(fin))
end

function main()
    target, pargs, figpath = parse_args(ARGS)
    cfg = load_config(ROOT)
    dirs = [ROOT]
    if isdir(target)
        pars = filter(f -> occursin(r"^yelmox.*\.nml$", f), readdir(target))
        isempty(pars) && error("no yelmox*.nml in $target")
        par = joinpath(target, pars[1])
        pushfirst!(dirs, target)
        figpath === nothing && (figpath = joinpath(target, "forcing.png"))
    else
        par = isfile(target) ? target : joinpath(ROOT, target)
        isfile(par) || error("par file $target not found")
        figpath === nothing && (figpath = "forcing-builder.png")
    end
    ov = parse_overrides(pargs)
    groups, order = resolve(cfg, par, ov)
    tl = timeline(groups)
    io = stdout

    println(io, "Par file: $par", isempty(ov) ? "" : "  (+ $(length(ov)) overrides)")
    println(io, "Domain: ", val(groups, "domain", "name", "?"), "; climate backend: ", val(groups, "comps", "climate", "?"), "; surface: ", val(groups, "comps", "surface_method", "?"))
    println(io)
    println(io, "Run control ([ctrl] run_step = $(tl.T))")
    @printf(io, "  time %s -> %s yr, dtt = %s, tstep_method = %s", fnum(tl.t0), fnum(tl.t1), fnum(val(groups, tl.T, "dtt", "?")), tl.method)
    tl.method == "const" && @printf(io, " (forcing held at t_f = %s)", fnum(tl.tc))
    tl.method == "cal" && @printf(io, " (calendar, t_f = t - %s)", fnum(tl.cref))
    println(io)
    println(io, "  restart = ", val(groups, "ctrl", "restart", "None"), "; init_method = ", val(groups, "sim", "init_method", "?"), "; relax = ", val(groups, "sim", "relax", "?"), ", opt = ", val(groups, "sim", "opt", "?"))
    if haskey(groups[tl.T].entries, "time_ref")
        for k in ("time_ref", "time_esm_ref", "time_hist", "time_proj")
            println(io, "  $k = ", join(string.(names_of(val(groups, tl.T, k, ""))), " - "))
        end
    end
    println(io)

    # indices
    idxnames = names_of(val(groups, "indices", "names", ""))
    println(io, "Indices: ", isempty(idxnames) ? "none" : join(idxnames, ", "))
    curves = Dict{String,Vector{Float64}}()
    for n in idxnames
        g = "idx_" * n
        haskey(groups, g) || (println(io, "  [$g]: missing"); continue)
        p = index_par(groups[g]); ser = load_index_series(p, dirs)
        f = index_series(p, tl.tf, tl.tf[1]; ser=ser)
        curves[n] = f
        describe_index(io, n, groups[g], p, ser, tl, f)
    end
    println(io)

    # snapshots / esm
    cl = string(val(groups, "comps", "climate", ""))
    if cl in ("snapshots", "rembo") && haskey(groups, "snap")
        snaps = names_of(val(groups, "snap", "snapshots", "")); axes = names_of(val(groups, "snap", "axes", ""))
        println(io, "Snapshots: ref = ", val(groups, "snap", "ref_name", "pd"), "; snapshots = ", join(snaps, ", "), "; f_p = ", fnum(val(groups, "snap", "f_p", 0.05)), ", lapse = ", join(fnum.(names_of(val(groups, "snap", "lapse", ""))), " "))
        for s in snaps
            g = "snap_" * s
            haskey(groups, g) || (println(io, "  [$g]: missing"); continue)
            bound = [k * " = " * join(names_of(val(groups, g, k, "")), " ") for k in ("tas", "pr", "zs", "to", "so") if !isempty(names_of(val(groups, g, k, "")))]
            println(io, "  [$g] ($(groups[g].source)): ", join(bound, "; "))
        end
        for a in axes
            g = "snap_axis_" * a
            haskey(groups, g) || (println(io, "  [$g]: missing"); continue)
            println(io, "  axis $a: lo = \"", val(groups, g, "lo", ""), "\", hi = \"", val(groups, g, "hi", ""), "\", index = ", val(groups, g, "index", ""), " (linear: X = X_ref + a (hi - lo); ratio: P = P_ref (1 + a (hi/lo - 1)))")
        end
        println(io)
    elseif cl == "esm"
        println(io, "ESM: experiment = ", val(groups, "esm", "experiment", "?"), ", esm_name = ", val(groups, "esm", "esm_name", "?"), ", use_esm = ", val(groups, "esm", "use_esm", "?"), ", use_hist = ", val(groups, "esm", "use_hist", "?"), ", use_proj = ", val(groups, "esm", "use_proj", "?"), ", par_file = ", val(groups, "esm", "par_file", "?"))
        println(io)
    end

    # modifiers
    println(io, "Climate modifiers")
    anoms = Dict{String,Vector{Float64}}()
    if val(groups, "anom", "active", false) == true
        ix = string(val(groups, "anom", "index", ""))
        fta = Float64(val(groups, "anom", "f_ta", 1.0)); fto = Float64(val(groups, "anom", "f_to", 0.0)); fso = Float64(val(groups, "anom", "f_so", 0.0))
        if cl == "rembo"
            @printf(io, "  anom: index %s, dT_sum = %g f, dT_ann = 1.3 dT_sum, dT_ocn = %g dT_ann\n", ix, fta, fto)
        else
            @printf(io, "  anom: index %s, dTa = %g f, dTo = %g f, dSo = %g f; pr -> pr exp(beta_p dTa)\n", ix, fta, fto, fso)
        end
        if haskey(curves, ix)
            anoms["dTa"] = curves[ix] .* fta
            fto != 0 && (anoms["dTo"] = (cl == "rembo" ? 1.3 .* curves[ix] .* fta .* fto : curves[ix] .* fto))
            fso != 0 && (anoms["dSo"] = curves[ix] .* fso)
            @printf(io, "    dTa over the run: %.3g at the start, %.3g at the end\n", anoms["dTa"][1], anoms["dTa"][end])
        else
            println(io, "    index \"$ix\" is not in [indices] names")
        end
    else
        println(io, "  anom: off")
    end
    val(groups, "ocn_fraction", "active", false) == true ?
        @printf(io, "  ocn_fraction: dTo = %g x mean(dT_sl) over %s\n", Float64(val(groups, "ocn_fraction", "f_to", 0.25)), val(groups, "ocn_fraction", "region", "all")) :
        println(io, "  ocn_fraction: off")
    val(groups, "taper", "active", false) == true ? println(io, "  taper: weight = ", val(groups, "taper", "weight", ""), ", outside = ", val(groups, "taper", "outside", "ref")) : println(io, "  taper: off")
    val(groups, "bias", "active", false) == true ? println(io, "  bias: dtas = ", val(groups, "bias", "dtas", ""), ", fpr = ", val(groups, "bias", "fpr", "")) : println(io, "  bias: off")
    println(io)

    # sea level, orbit, ice_bnd
    slm = string(val(groups, "barysealevel", "method", "fastiso"))
    sl = nothing
    if slm == "file"
        path = string(val(groups, "barysealevel", "sl_path", "")); vn = string(val(groups, "barysealevel", "sl_name", "none"))
        for d in dirs
            isfile(joinpath(d, path)) || continue
            ser = load_series(joinpath(d, path); var=vn)
            sl = [interp(ser, x) for x in tl.tf]
            break
        end
        println(io, "Sea level: from $path", sl === nothing ? " (file not found)" : @sprintf(" (%.4g m at the start, %.4g m at the end)", sl[1], sl[end]))
    else
        println(io, "Sea level: $slm", slm == "const" ? @sprintf(" (%s m)", fnum(val(groups, "barysealevel", "bsl_init", 0.0))) : " (from the FastIsostasy ice volume)")
    end
    om = string(val(groups, "orbit", "method", "time"))
    println(io, "Orbit: ", om == "clock" ? "at the forcing clock" : om == "time" ? "fixed at $(fnum(val(groups, "orbit", "time", 0.0))) yr" : "fixed elements ecc = $(val(groups, "orbit", "ecc", "")), obl = $(val(groups, "orbit", "obl", "")), per = $(val(groups, "orbit", "per", ""))")
    if val(groups, "ice_bnd", "active", false) == true
        println(io, "Prescribed ice: ", val(groups, "ice_bnd", "method", "file"), ", H_ice = ", val(groups, "ice_bnd", "H_ice", ""), val(groups, "ice_bnd", "method", "") == "scaled" ? ", index $(val(groups, "ice_bnd", "index", "")) x H_ice($(fnum(val(groups, "ice_bnd", "slice_time", 0.0))))" : "")
    end
    println(io)

    # the sources of the forcing values
    println(io, "Sources: ", join(["$(g) ($(groups[g].source == "new" ? "defaults only" : groups[g].source))" for g in order if groups[g].source != "base"], ", "))
    overridden = [g * "." * k for g in order for k in groups[g].order if groups[g].entries[k].source == "override"]
    isempty(overridden) || println(io, "Overrides: ", join(overridden, ", "))

    # figure
    (figpath == "none") && return
    CM = Base.require(Main, :CairoMakie)
    Base.invokelatest() do
        fig = CM.Figure(size=(800, 220 * (1 + length(curves) + (isempty(anoms) ? 0 : 1) + (sl === nothing ? 0 : 1))))
        row = 1
        ax = CM.Axis(fig[row, 1], title="Forcing clock ($(tl.method))", xlabel="t [yr]", ylabel="t_f [yr]"); CM.lines!(ax, tl.t, tl.tf); row += 1
        for n in idxnames
            haskey(curves, n) || continue
            ax = CM.Axis(fig[row, 1], title="Index $n", xlabel="t [yr]", ylabel=n); CM.lines!(ax, tl.t, curves[n]); row += 1
        end
        if !isempty(anoms)
            ax = CM.Axis(fig[row, 1], title="Homogeneous anomalies [anom]", xlabel="t [yr]", ylabel="anomaly")
            for (k, v) in sort(collect(anoms)); CM.lines!(ax, tl.t, v, label=k); end
            CM.axislegend(ax, position=:lt); row += 1
        end
        if sl !== nothing
            ax = CM.Axis(fig[row, 1], title="Sea level [barysealevel]", xlabel="t [yr]", ylabel="z_sl [m]"); CM.lines!(ax, tl.t, sl); row += 1
        end
        CM.save(figpath, fig)
    end
    println(io, "Figure: $figpath")
end

main()

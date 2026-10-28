# check_forcing.jl -- the forcing of a yelmox_bipolar run, checked from its
# output against its own par file (the groups not in it: the library
# input/forcing_db.nml). Run from the yelmox root:
#
#     julia --project=yelmox_bipolar yelmox_bipolar/check_forcing.jl <run> [<reference run>]
#
# 1. gmt (forcing_ts.nc) follows [coef_gmt] (const or ramp-time);
# 2. each scaled coefficient is factor x its source (attributes kind, coef,
#    factor);
# 3. per active domain (its folder's snap.nc): term_<name> is the value of its
#    coefficient; relative to the first record, tsl_ann moves by the uniform
#    tas terms (f_ta x coefficient), pr_ann by exp(f_p x the uniform pr
#    terms), to_ann by the uniform to terms, or is the replace term's value;
#    a derived air_anomaly coefficient moves by factor x the mean tsl_ann
#    anomaly;
# 4. the box model (<obm_name>.nc, atm2obm): theta_i = theta_i_init + its
#    coefficient, phin/phit linear in the global coefficient; obm_north/
#    obm_south are the box temperature in K (minus its initial value with
#    anomaly = True);
# 5. with a reference run (a restart of it): every value equal at the common
#    times.
# Prints each check; exits 1 if one fails. Written for the forcing tests of
# run_tests.sh (ice sheets frozen: the elevation, and so the lapse-rate
# correction, constant).

using NCDatasets, Printf
include(joinpath(@__DIR__, "..", "scripts", "forcing-builder", "Namelist.jl"))
using .Namelist

const NFAIL = Ref(0)

function check(ok::Bool, what::AbstractString)
    println(ok ? "  ok    " : "  FAIL  ", what)
    ok || (NFAIL[] += 1)
    return ok
end

# ----- the run's parameters ----------------------------------------------

"The run's par file (the one *.nml of the run folder) and the library."
function load_par(run)
    nmls = filter(f -> endswith(f, ".nml"), readdir(run))
    length(nmls) == 1 || error("$run: expected one par file (*.nml), found $(nmls)")
    files = [group_dict(read_namelist(joinpath(run, nmls[1])))]
    lib = joinpath(run, "input", "forcing_db.nml")
    isfile(lib) && push!(files, group_dict(read_namelist(lib)))
    return files
end

"A key of a group, from the first file that has the group (as the model reads
it); `default` when the group or the key is absent."
function par(P, group, key, default=nothing)
    for d in P
        haskey(d, group) || continue
        for k in d[group].keys
            k.name == key && return parse_value(k.raw)
        end
        break
    end
    default === nothing && error("[$group] $key not found")
    return default
end

aslist(v) = v isa AbstractVector ? string.(v) : (v == "" ? String[] : [string(v)])

# ----- output ------------------------------------------------------------

"A variable as a Float64 array of its dimensions, fill values as NaN."
readvar(ds, v) = Float64.(coalesce.(Array(ds[v]), NaN))

"The record indices of the times `t` in `tf` (all must be there)."
records(tf, t) = [findfirst(==(x), tf) for x in t]

maxabs(x) = isempty(x) ? 0.0 : maximum(abs, x)

# ----- checks ------------------------------------------------------------

function check_gmt(P, ts, t0)
    method = par(P, "coef_gmt", "method")
    gmt = readvar(ts, "gmt")
    t = readvar(ts, "time")
    f_min = par(P, "coef_gmt", "f_min", 0.0)
    f_max = par(P, "coef_gmt", "f_max", 1.0)
    if method == "const"
        expect = fill(Float64(f_min), length(t))
    elseif method == "ramp-time"
        df_sign = par(P, "coef_gmt", "df_sign", 1.0)
        dt_ramp = par(P, "coef_gmt", "dt_ramp", 1000.0)
        f_start = df_sign > 0 ? f_min : f_max
        rate = abs(f_max - f_min) / dt_ramp
        expect = clamp.(f_start .+ df_sign * rate .* (t .- t0), f_min, f_max)
    else
        println("  --    gmt: method $method not checked")
        return
    end
    check(maxabs(gmt .- expect) < 1e-6,
          @sprintf("gmt follows [coef_gmt] %s (max |diff| %.1e K)", method, maxabs(gmt .- expect)))
end

function check_scaled(ts)
    for v in keys(ts)
        get(ts[v].attrib, "kind", "") == "scaled" || continue
        src, factor = ts[v].attrib["coef"], ts[v].attrib["factor"]
        d = maxabs(readvar(ts, v) .- factor .* readvar(ts, src))
        check(d < 1e-6, @sprintf("%s = %.2f x %s (max |diff| %.1e)", v, factor, src, d))
    end
end

"The terms of a domain: name, coefficient, pattern, fields, f_ta, f_to."
function domain_terms(P, sfx)
    terms = []
    for name in aslist(par(P, "forcing$sfx", "terms", ""))
        g = "term_$name"
        push!(terms, (name = name, coef = string(par(P, g, "coef", "")),
                      pattern = string(par(P, g, "pattern", "snapshots")),
                      fields = aslist(par(P, g, "fields")),
                      f_ta = Float64(par(P, g, "f_ta", 1.0)), f_to = Float64(par(P, g, "f_to", 1.0))))
    end
    return terms
end

function check_domain(P, run, sfx, ts)
    name = par(P, "domain$sfx", "name")
    fnm = joinpath(run, name, "snap.nc")
    println("$name ($fnm)")
    isfile(fnm) || return check(false, "$fnm exists")
    terms = domain_terms(P, sfx)
    sn = NCDataset(fnm)
    tf, t = readvar(ts, "time"), readvar(sn, "time")
    kf = records(tf, t)
    any(isnothing, kf) && return check(false, "the times of snap.nc are in forcing_ts.nc")
    coef(c) = readvar(ts, c)[kf]

    # Each term reads its coefficient.
    for tm in terms
        tm.coef == "" && continue
        d = maxabs(readvar(sn, "term_$(tm.name)") .- coef(tm.coef))
        check(d < 1e-6, @sprintf("term_%s = %s (max |diff| %.1e)", tm.name, tm.coef, d))
    end

    # The atmosphere: the uniform terms on tas and pr (another pattern on
    # them: not checked).
    sum_uniform(field, f) = sum((f(tm) .* coef(tm.coef) for tm in terms
                                 if field in tm.fields && tm.pattern == "uniform"), init = zeros(length(t)))
    others(field) = any(field in tm.fields && tm.pattern != "uniform" for tm in terms)
    tsl = readvar(sn, "tsl_ann")
    dtsl = tsl .- tsl[:, :, 1:1]
    if !others("tas")
        s = sum_uniform("tas", tm -> tm.f_ta)
        d = maxabs(filter(isfinite, dtsl .- reshape(s .- s[1], 1, 1, :)))
        check(d < 2e-3, @sprintf("tsl_ann moves by the uniform tas terms (max |diff| %.1e K, tsl_ann written to 1e-3 K)", d))
    end
    if !others("pr") && haskey(sn, "pr_ann")
        f_p = Float64(par(P, "snapshots$sfx", "f_p"))
        s = sum_uniform("pr", tm -> tm.f_ta)
        pr = readvar(sn, "pr_ann")
        ratio = pr ./ pr[:, :, 1:1]
        ok = (pr[:, :, 1:1] .> 1e-3) .& isfinite.(ratio)
        expect = reshape(exp.(f_p .* (s .- s[1])), 1, 1, :)
        d = maxabs((ratio .- expect)[repeat(ok, 1, 1, length(t))])
        check(d < 2e-3, @sprintf("pr_ann scales by exp(f_p = %.3f x the uniform pr terms) (max |rel diff| %.1e)", f_p, d))
    end

    # A derived air_anomaly coefficient read by this domain: factor x the
    # mean tsl_ann anomaly of the climate grid.
    for tm in terms
        tm.coef == "" && continue
        get(ts[tm.coef].attrib, "from", "") == "air_anomaly" || continue
        factor = Float64(par(P, "coef_$(tm.coef)", "factor", 1.0))
        c = coef(tm.coef)
        mean_dtsl = [sum(filter(isfinite, dtsl[:, :, k])) / count(isfinite, dtsl[:, :, k]) for k in eachindex(t)]
        d = maxabs((c .- c[1]) .- factor .* mean_dtsl)
        check(d < 2e-3, @sprintf("%s moves by %.2f x the mean tsl_ann anomaly (max |diff| %.1e K)", tm.coef, factor, d))
    end

    # The ocean: a replace term sets to_ann, the uniform terms add to it; on
    # the cells with an ocean temperature at the first record (no fill value).
    if any("to" in tm.fields for tm in terms) && haskey(sn, "to_ann")
        to = readvar(sn, "to_ann")
        wet = repeat((to[:, :, :, 1:1] .> 250) .& (to[:, :, :, 1:1] .< 320), 1, 1, 1, length(t))
        u = sum_uniform("to", tm -> tm.f_to)
        rep = [tm for tm in terms if "to" in tm.fields && tm.pattern == "replace"]
        if length(rep) == 1
            expect = coef(rep[1].coef) .+ u
            d = maxabs((to .- reshape(expect, 1, 1, 1, :))[wet])
            check(d < 2e-3, @sprintf("to_ann = %s (replace) + the uniform to terms (max |diff| %.1e K, %d wet cells)",
                                     rep[1].coef, d, count(wet[:, :, :, 1])))
        elseif isempty(rep)
            d = maxabs((to .- to[:, :, :, 1:1] .- reshape(u .- u[1], 1, 1, 1, :))[wet])
            check(d < 2e-3, @sprintf("to_ann moves by the uniform to terms (max |diff| %.1e K, %d wet cells)",
                                     d, count(wet[:, :, :, 1])))
        end
    end
    close(sn)
end

function check_obm(P, run, ts)
    obm_name = par(P, "ctrl", "obm_name")
    fnm = joinpath(run, "$obm_name.nc")
    println("box model ($fnm)")
    isfile(fnm) || return check(false, "$fnm exists")
    ob = NCDataset(fnm)
    tf, t = readvar(ts, "time"), readvar(ob, "time")
    kf = records(tf, t)
    any(isnothing, kf) && return check(false, "the times of $obm_name.nc are in forcing_ts.nc")
    coef(c) = readvar(ts, c)[kf]
    np(key) = Float64(par(P, obm_name, key))

    if par(P, "ctrl", "atm2obm", false) == true
        for (box, key) in (("n", "coef_north"), ("t", "coef_tropics"), ("s", "coef_south"))
            c = par(P, "obm_forcing", key)
            d = maxabs(readvar(ob, "theta$box") .- (np("theta$(box)_init") .+ coef(c)))
            check(d < 1e-5, @sprintf("theta%s = theta%s_init + %s (max |diff| %.1e K)", box, box, c, d))
        end
        g = coef(par(P, "obm_forcing", "coef"))
        d = maxabs(readvar(ob, "phin") .- (np("phin_init") .+ np("nh_hydro_sensitivity") * np("nh_temp_constant") .* g))
        check(d < 1e-7, @sprintf("phin linear in the global coefficient (max |diff| %.1e Sv)", d))
        d = maxabs(readvar(ob, "phit") .- (np("phit_init") .+ np("sh_hydro_sensitivity") * np("sh_temp_constant") .* g))
        check(d < 1e-7, @sprintf("phit linear in the global coefficient (max |diff| %.1e Sv)", d))
    end

    for (c, box) in (("obm_north", "tn"), ("obm_south", "ts"))
        haskey(ts, c) || continue
        anomaly = par(P, "coef_$c", "anomaly", true) == true
        T = readvar(ob, box) .+ 273.15
        expect = anomaly ? T .- (np("$(box)_init") + 273.15) : T
        d = maxabs(coef(c) .- expect)
        check(d < 1e-4, @sprintf("%s = %s%s in K (max |diff| %.1e K)", c, box, anomaly ? " - $(box)_init" : "", d))
    end
    close(ob)
end

"Every variable of a file equal to the reference's at the common times."
function check_same(fnm, fref; vars = nothing)
    (isfile(fnm) && isfile(fref)) || return
    a, b = NCDataset(fnm), NCDataset(fref)
    ta, tb = readvar(a, "time"), readvar(b, "time")
    common = intersect(ta, tb)
    ka, kb = records(ta, common), records(tb, common)
    worst, vworst = 0.0, ""
    for v in (vars === nothing ? keys(a) : vars)
        (v == "time" || !haskey(b, v) || !("time" in dimnames(a[v]))) && continue
        x, y = readvar(a, v), readvar(b, v)
        nd = ndims(x)
        xa = selectdim(x, nd, ka); yb = selectdim(y, nd, kb)
        d = maxabs(filter(isfinite, xa .- yb))
        d > worst && ((worst, vworst) = (d, v))
    end
    check(worst < 1e-6, @sprintf("%s = reference at %d common times (max |diff| %.1e%s)",
                                 basename(fnm), length(common), worst, vworst == "" ? "" : " in $vworst"))
    close(a); close(b)
end

function print_table(ts, names)
    names = [n for n in names if haskey(ts, n)]
    println("  ", rpad("time", 8), join(lpad.(names, 11)))
    cols = [readvar(ts, n) for n in names]
    for (k, t) in enumerate(readvar(ts, "time"))
        println("  ", rpad(@sprintf("%.0f", t), 8), join(@sprintf("%11.4f", c[k]) for c in cols))
    end
end

# ----- main --------------------------------------------------------------

function main(args)
    (1 <= length(args) <= 2) || error("usage: check_forcing.jl <run> [<reference run>]")
    run = args[1]
    ref = length(args) == 2 ? args[2] : nothing
    P = load_par(run)

    println("== $run", ref === nothing ? "" : "  (reference $ref)")
    ts = NCDataset(joinpath(run, "forcing_ts.nc"))
    print_table(ts, ["gmt", "dTa_north", "dTa_south", "dTo", "obm_north", "obm_south", "bsl"])

    println("forcing_ts.nc")
    # The ramp starts at the start of the run, or of its reference (a restart
    # continues the ramp).
    t0 = ref === nothing ? readvar(ts, "time")[1] :
         NCDataset(f -> readvar(f, "time")[1], joinpath(ref, "forcing_ts.nc"))
    check_gmt(P, ts, t0)
    check_scaled(ts)

    for sfx in ("_north", "_south")
        par(P, "ctrl", "active$sfx") == true && check_domain(P, run, sfx, ts)
    end
    par(P, "ctrl", "active_obm") == true && check_obm(P, run, ts)

    if ref !== nothing
        println("against the reference")
        check_same(joinpath(run, "forcing_ts.nc"), joinpath(ref, "forcing_ts.nc"))
        obm_name = par(P, "ctrl", "obm_name")
        check_same(joinpath(run, "$obm_name.nc"), joinpath(ref, "$obm_name.nc"))
        for sfx in ("_north", "_south")
            par(P, "ctrl", "active$sfx") == true || continue
            name = par(P, "domain$sfx", "name")
            check_same(joinpath(run, name, "snap.nc"), joinpath(ref, name, "snap.nc"),
                       vars = ["tsl_ann", "pr_ann", "to_ann"])
        end
    end
    close(ts)

    println(NFAIL[] == 0 ? "check_forcing: PASS" : "check_forcing: FAIL ($(NFAIL[]) checks)")
    exit(NFAIL[] == 0 ? 0 : 1)
end

main(ARGS)

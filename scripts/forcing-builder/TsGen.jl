# TsGen.jl — the methods of fesm-utils tsgen (src/tsgen.f90) as pure
# functions, to tabulate a coefficient [coef_<name>] on the forcing time without
# running the model (a derived coefficient is set by the model: not here). Same
# formulas as docs/forcing-builder (app.js).
#
# Open-loop methods (eval_open), with tau = max(0, t_elapsed - dt_init) and
# f_start = f_min (df_sign > 0) or f_max:
#   const, ramp-slope : f = clamp(f_start + df_sign*df_dt_max*tau, f_min, f_max)
#   ramp-time         : f = clamp(f_start + df_sign*|f_max - f_min|/dt_ramp*tau)
#   ramp-time-step    : f_start -> f_ext over dt_ramp -> f_conv over dt_conv -> hold
#   sin               : f = (f_max + f_min)/2 + (f_max - f_min)/2 sin(2 pi t_elapsed/dt_ramp)
#   series            : f = S(t) (normalised: (S - S0)/(S1 - S0))
#   feedback (exp, PI42, H312b, H312PID, H321PID, PID1): no analytic form, f_start
# then the window modifier (apply_window) and noise (sigma).
module TsGen

export Series, load_series, interp, series_mean, normalize, eval_open,
       apply_window, tabulate, is_feedback, is_open_loop

struct Series
    time::Vector{Float64}
    f::Vector{Float64}     # the first channel (the channel mean for several)
    header::String
    nc::Int                # number of channels in the file
end

const FEEDBACK = ("exp", "PI42", "H312b", "H312PID", "H321PID", "PID1")
const OPEN = ("const", "ramp-slope", "ramp-time", "ramp-time-step", "sin")

is_feedback(method) = method in FEEDBACK
is_open_loop(method) = method in OPEN

"Linear interpolation, held at the ends (interp_linear of fesm-utils)."
function interp(ser::Series, t::Real)
    n = length(ser.time)
    n == 0 && return NaN
    t <= ser.time[1] && return ser.f[1]
    t >= ser.time[end] && return ser.f[end]
    i = searchsortedlast(ser.time, t)
    t0, t1 = ser.time[i], ser.time[i+1]
    w = (t - t0) / (t1 - t0)
    return ser.f[i] * (1 - w) + ser.f[i+1] * w
end

"series_mean: the value at one time, or the mean of the samples in a window."
function series_mean(ser::Series, tw)
    tw = collect(Float64, tw)
    length(tw) == 1 && return interp(ser, tw[1])
    lo, hi = minimum(tw), maximum(tw)
    sel = [f for (t, f) in zip(ser.time, ser.f) if lo <= t <= hi]
    isempty(sel) && return interp(ser, (lo + hi) / 2)
    return sum(sel) / length(sel)
end

"The series rescaled to 0 at norm_t0 and 1 at norm_t1 (one time or a window each)."
function normalize(ser::Series, t0, t1)
    v0 = series_mean(ser, t0)
    v1 = series_mean(ser, t1)
    return Series(ser.time, (ser.f .- v0) ./ (v1 - v0), ser.header, ser.nc)
end

function load_series_ascii(path::AbstractString)
    time = Float64[]; f = Float64[]; header = String[]; nc = 0
    for line in eachline(path)
        s = strip(line)
        isempty(s) && continue
        if startswith(s, "#") || startswith(s, "!")
            push!(header, String(strip(s[2:end])))
            continue
        end
        toks = split(s)
        vals = [tryparse(Float64, t) for t in toks]
        if any(v -> v === nothing, vals)
            push!(header, String(s))   # a text header line such as "time alpha"
            continue
        end
        length(vals) < 2 && continue
        nc = max(nc, length(vals) - 1)
        push!(time, vals[1]); push!(f, vals[2])
    end
    return Series(time, f, join(header, " "), nc)
end

"""
Load a tsgen series: an ASCII table (time, f, ...; lines starting with # or !
skipped) or a netCDF variable `var` on `time_name`.
"""
function load_series(path::AbstractString; var::AbstractString="", time_name::AbstractString="time")
    ext = lowercase(splitext(path)[2])
    if ext in (".nc", ".nc4", ".cdf", ".netcdf")
        NC = Base.require(Main, :NCDatasets)
        ds = NC.Dataset(path)
        try
            t = Float64.(collect(ds[time_name][:]))
            v = ds[var]
            data = v[:]
            vals = ndims(data) == 1 ? Float64.(collect(data)) :
                   vec(sum(Float64.(collect(data)), dims=2)) ./ size(data, 2)
            return Series(t, vals, "$(path):$(var)", ndims(data) == 1 ? 1 : size(data, 2))
        finally
            close(ds)
        end
    end
    return load_series_ascii(path)
end

getpar(par, key, default) = haskey(par, key) ? par[key] : default

"The open-loop value at the elapsed time since the series was initialised."
function eval_open(par::AbstractDict, t_elapsed::Real)
    method = getpar(par, "method", "const")
    dt_init = getpar(par, "dt_init", 0.0)
    f_min = getpar(par, "f_min", 0.0); f_max = getpar(par, "f_max", 1.0)
    sgn = getpar(par, "df_sign", 1.0)
    tau = max(0.0, t_elapsed - dt_init)
    f_start = sgn > 0 ? f_min : f_max
    clamp01(x) = min(max(x, f_min), f_max)
    if method in ("const", "ramp-slope")
        return clamp01(f_start + sgn * getpar(par, "df_dt_max", 0.0) * tau)
    elseif method == "ramp-time"
        rate = abs(f_max - f_min) / getpar(par, "dt_ramp", 1000.0)
        return clamp01(f_start + sgn * rate * tau)
    elseif method == "ramp-time-step"
        dt_ramp = getpar(par, "dt_ramp", 1000.0); dt_conv = getpar(par, "dt_conv", 1000.0)
        f_conv = getpar(par, "f_conv", 0.0)
        f_ext = sgn > 0 ? f_max : f_min
        if tau <= dt_ramp
            return f_start + (f_ext - f_start) * (tau / dt_ramp)
        elseif tau <= dt_ramp + dt_conv
            return f_ext + (f_conv - f_ext) * ((tau - dt_ramp) / dt_conv)
        else
            return f_conv
        end
    elseif method == "sin"
        p = getpar(par, "dt_ramp", 1000.0)
        te = t_elapsed <= dt_init ? 0.0 : t_elapsed
        return 0.5 * (f_max + f_min) + 0.5 * (f_max - f_min) * sin(2pi * te / p)
    else
        return f_start      # feedback methods: no analytic form
    end
end

"The window modifier: fw = f + w(t)*(win_fac - 1)*(f - win_ref)."
function apply_window(par::AbstractDict, t::Real, f::Real)
    t0 = getpar(par, "win_t0", 0.0); t1 = getpar(par, "win_t1", 0.0)
    t1 <= t0 && return f
    (t <= t0 || t >= t1) && return f
    ref = getpar(par, "win_ref", 0.0)
    getpar(par, "win_above", false) && f <= ref && return f
    taper = getpar(par, "win_taper", 0.0)
    w = 1.0
    taper > 0 && (w = min(1.0, (t - t0) / taper, (t1 - t) / taper))
    return f + w * (getpar(par, "win_fac", 1.0) - 1.0) * (f - ref)
end

"""
The coefficient on the forcing time: `times` the forcing times of the run and
`t_init` the forcing time at initialisation (the elapsed time of the open-loop
methods counts from it). `ser` is the loaded (normalised) series for "series".
Returns the mean value (no noise), windowed.
"""
function tabulate(par::AbstractDict, times::AbstractVector, t_init::Real; ser::Union{Series,Nothing}=nothing)
    method = getpar(par, "method", "const")
    out = similar(collect(Float64, times))
    for (i, t) in enumerate(times)
        f = method == "series" ? (ser === nothing ? NaN : interp(ser, t)) : eval_open(par, t - t_init)
        out[i] = apply_window(par, t, f)
    end
    return out
end

end # module

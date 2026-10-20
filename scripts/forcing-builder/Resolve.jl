# Resolve.jl — the configuration a run sees: a base par file, the named groups
# it references from the library (.runme/info.json par_library), the schema
# defaults (par_defaults: input/yelmox_defaults.nml) and the runme -p overrides,
# with the source of every value (base, library, default, override).
module Resolve

using ..Namelist
using JSON

export Config, load_config, resolve, parse_overrides, template, RGroup, Entry, value_of

const TEMPLATES = [
    (r"^idx_\w+$", "idx"),
    (r"^snap_field_\w+$", "snap_field"),
    (r"^snap_axis_\w+$", "snap_axis"),
    (r"^snap_\w+$", "snap_snapshot"),
]

struct Entry
    raw::String
    value
    source::String      # base | library | default | override
    comment::String
    pre::String
end

struct RGroup
    name::String
    source::String      # base | library | new
    template::Union{String,Nothing}
    order::Vector{String}
    entries::Dict{String,Entry}
end

struct Config
    root::String
    defaults::Dict{String,Group}
    library::Dict{String,Group}
    library_files::Vector{String}
    defaults_files::Vector{String}
end

"The defaults group of a named group (its own name, or the template of its kind)."
function template(name::AbstractString, defaults::Dict{String,Group})
    haskey(defaults, name) && return String(name)
    for (re, t) in TEMPLATES
        occursin(re, name) && return t
    end
    return nothing
end

"Read .runme/info.json (par_defaults, par_library) and the files it names."
function load_config(root::AbstractString)
    info = JSON.parsefile(joinpath(root, ".runme/info.json"))
    dfiles = String.(get(info, "par_defaults", String[]))
    lfiles = String.(get(info, "par_library", String[]))
    defaults = Dict{String,Group}()
    for f in dfiles
        endswith(f, "yelmox_defaults.nml") || continue   # the yelmo groups are not forcing
        merge!(defaults, group_dict(read_namelist(joinpath(root, f))))
    end
    library = Dict{String,Group}()
    for f in lfiles
        for g in read_namelist(joinpath(root, f))
            haskey(library, g.name) || (library[g.name] = g)
        end
    end
    return Config(String(root), defaults, library, lfiles, dfiles)
end

"""
The -p tokens of runme as (key, tokens, vector): `g.k=v` a scalar, `g.k==a,b`
a vector. Ensemble forms (comma lists after a single =, ranges) are rejected.
"""
function parse_overrides(args::AbstractVector{<:AbstractString})
    out = Tuple{String,Vector{String},Bool}[]
    for a in args
        if occursin("==", a)
            k, v = split(a, "==", limit=2)
            push!(out, (String(k), String.(split(v, ",")), true))
        else
            occursin("=", a) || error("override '$a': expected key=value")
            k, v = split(a, "=", limit=2)
            occursin(",", v) && error("override '$a': an ensemble list; use == for a vector")
            push!(out, (String(k), [String(v)], false))
        end
    end
    return out
end

"Coerce the text of an override to the type of the reference value."
function coerce(tok::AbstractString, ref)
    ref isa Bool && return lowercase(tok) in ("true", ".true.", "t", "1")
    ref isa Integer && (x = tryparse(Int, tok); x !== nothing && return x)
    ref isa Real && (x = tryparse(Float64, replace(lowercase(tok), "d" => "e")); x !== nothing && return x)
    return String(strip(tok, ['"', '\'']))
end

function override_value(toks, vector::Bool, ref)
    refel = ref isa AbstractVector ? (isempty(ref) ? "" : ref[1]) : ref
    vals = [coerce(t, refel) for t in toks]
    (vector || ref isa AbstractVector) && return vals
    return vals[1]
end

value_of(rg::RGroup, key::AbstractString, default=nothing) =
    haskey(rg.entries, key) ? rg.entries[key].value : default

"Names listed in a vector-of-strings key (\"\" entries dropped)."
function names_of(rg::Union{RGroup,Nothing}, key::AbstractString)
    rg === nothing && return String[]
    v = value_of(rg, key, "")
    v isa AbstractVector && return [string(x) for x in v if !isempty(string(x))]
    s = string(v)
    return isempty(s) ? String[] : [s]
end

"""
Resolve `base_path` with `overrides` (from parse_overrides). Returns the groups
(Dict name => RGroup) and their order: the base's groups, then the named groups
the configuration references (indices, snapshots, axes, fields), from the base,
the library, or new from the template.
"""
function resolve(cfg::Config, base_path::AbstractString, overrides=Tuple{String,Vector{String},Bool}[])
    base_groups = read_namelist(base_path)
    base = group_dict(base_groups)
    out = Dict{String,RGroup}()
    order = String[]

    function group!(name::AbstractString)
        haskey(out, name) && return out[name]
        tname = template(name, cfg.defaults)
        if haskey(base, name)
            g = base[name]; src = "base"
        elseif haskey(cfg.library, name)
            g = cfg.library[name]; src = "library"
        elseif tname !== nothing
            g = Group(String(name), "", Key[]); src = "new"
        else
            error("group [$name]: not in the par file, the library or the defaults")
        end
        entries = Dict{String,Entry}()
        ord = String[]
        if tname !== nothing
            for k in cfg.defaults[tname].keys
                entries[k.name] = Entry(k.raw, parse_value(k.raw), "default", k.comment, k.pre)
                push!(ord, k.name)
            end
        end
        for k in g.keys
            if haskey(entries, k.name)
                c = isempty(k.comment) ? entries[k.name].comment : k.comment
                p = isempty(k.pre) ? entries[k.name].pre : k.pre
            else
                push!(ord, k.name); c = k.comment; p = k.pre
            end
            entries[k.name] = Entry(k.raw, parse_value(k.raw), src, c, p)
        end
        rg = RGroup(String(name), src, tname, ord, entries)
        out[name] = rg
        push!(order, name)
        return rg
    end

    for g in base_groups
        group!(g.name)
    end
    for (key, toks, vector) in overrides
        occursin(".", key) || error("override '$key': expected group.key")
        gname, kname = split(key, ".", limit=2)
        rg = group!(gname)
        haskey(rg.entries, kname) || error("override '$key': key not in [$gname] nor its defaults")
        e = rg.entries[kname]
        v = override_value(toks, vector, e.value)
        rg.entries[kname] = Entry(format_value(v), v, "override", e.comment, e.pre)
    end
    for n in names_of(get(out, "indices", nothing), "names")
        group!("idx_" * n)
    end
    if haskey(out, "snap")
        for n in names_of(out["snap"], "snapshots"); group!("snap_" * n); end
        for n in names_of(out["snap"], "axes");      group!("snap_axis_" * n); end
        for n in names_of(out["snap"], "fields");    group!("snap_field_" * n); end
    end
    return out, order
end

end # module

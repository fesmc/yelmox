# Namelist.jl — reader of the Fortran namelist files of yelmox, as fesm-utils
# nml.f90 reads them (case-sensitive names, one key per line, an inline comment
# after a ! outside quotes), keeping the comments: the schema of the forcing
# groups is input/yelmox_defaults.nml with its comments.
#
#   groups = read_namelist(path)         # Vector{Group}, in file order
#   parse_value("\"pd\" \"lgm\"")         # ["pd", "lgm"]
#   format_value(["pd", "lgm"])          # "\"pd\" \"lgm\""
module Namelist

export Group, Key, read_namelist, parse_value, format_value, group_dict

struct Key
    name::String
    raw::String        # the value as written, without the comment
    comment::String    # the inline comment (and its continuation lines)
    pre::String        # comment lines before the key inside the group
end

struct Group
    name::String
    comment::String    # comment lines right after &name, before the first key
    keys::Vector{Key}
end

"The value without its inline comment (a ! outside quotes), and the comment."
function split_comment(s::AbstractString)
    q = nothing
    for (i, c) in enumerate(s)
        if q !== nothing
            c == q && (q = nothing)
        elseif c == '"' || c == '\''
            q = c
        elseif c == '!'
            return (strip(s[1:prevind(s, i)]), strip(s[nextind(s, i):end]))
        end
    end
    return (strip(s), "")
end

"Column of the ! that starts a comment (1-based), or 0."
function comment_column(line::AbstractString)
    q = nothing
    for (i, c) in enumerate(line)
        if q !== nothing
            c == q && (q = nothing)
        elseif c == '"' || c == '\''
            q = c
        elseif c == '!'
            return i
        end
    end
    return 0
end

function read_namelist(path::AbstractString)
    groups = Group[]
    name = nothing
    gcomment = String[]
    keys = Key[]
    pre = String[]
    lastcol = 0          # column of the inline ! of the previous key
    for line in eachline(path)
        s = strip(line)
        if isempty(s)
            lastcol = 0
            continue
        end
        if startswith(s, "&")
            name = String(strip(s[2:end]))
            gcomment = String[]; keys = Key[]; pre = String[]; lastcol = 0
            continue
        end
        if startswith(s, "/")
            name === nothing && continue
            push!(groups, Group(name, join(gcomment, "\n"), keys))
            name = nothing
            continue
        end
        name === nothing && continue
        if startswith(s, "!")
            text = String(strip(s[2:end]))
            col = comment_column(line)
            if !isempty(keys) && lastcol > 0 && col == lastcol
                # continuation of the previous key's inline comment
                k = keys[end]
                keys[end] = Key(k.name, k.raw, k.comment * "\n" * text, k.pre)
            elseif isempty(keys) && isempty(pre)
                push!(gcomment, text)
            else
                push!(pre, text)
                lastcol = 0
            end
            continue
        end
        occursin("=", s) || continue
        k, v = split(s, "=", limit=2)
        raw, comment = split_comment(v)
        raw = String(strip(rstrip(raw, ',')))
        push!(keys, Key(String(strip(k)), raw, comment, join(pre, "\n")))
        pre = String[]
        lastcol = isempty(comment) ? 0 : comment_column(line)
    end
    return groups
end

"Groups as a Dict name => Group (the first copy wins)."
function group_dict(groups::Vector{Group})
    d = Dict{String,Group}()
    for g in groups
        haskey(d, g.name) || (d[g.name] = g)
    end
    return d
end

const RE_TOKEN = r"\"([^\"]*)\"|'([^']*)'|([^\s,\"']+)"

"Tokens of a value: (text, quoted)."
function tokens(raw::AbstractString)
    out = Tuple{String,Bool}[]
    for m in eachmatch(RE_TOKEN, raw)
        if m.captures[1] !== nothing
            push!(out, (m.captures[1], true))
        elseif m.captures[2] !== nothing
            push!(out, (m.captures[2], true))
        else
            push!(out, (m.captures[3], false))
        end
    end
    return out
end

function parse_scalar(tok::AbstractString, quoted::Bool)
    quoted && return String(tok)
    t = lowercase(tok)
    (t in ("true", ".true.", "t")) && return true
    (t in ("false", ".false.", "f")) && return false
    tn = replace(t, "d" => "e")
    occursin(r"^[-+]?\d+$", tn) && return parse(Int, tn)
    x = tryparse(Float64, tn)
    x === nothing ? String(tok) : x
end

"""
The Julia value of a namelist value: a String, Bool, Int or Float64, or a
Vector of them when the value has several tokens. `""` is the empty string.
"""
function parse_value(raw::AbstractString)
    toks = tokens(raw)
    isempty(toks) && return ""
    vals = [parse_scalar(t, q) for (t, q) in toks]
    length(vals) == 1 && return vals[1]
    any(last, toks) && return [string(v) for v in vals]
    return vals
end

format_scalar(v::Bool) = v ? "True" : "False"
format_scalar(v::Integer) = string(v)
format_scalar(v::AbstractFloat) = string(v)
format_scalar(v::AbstractString) = "\"" * v * "\""

"The namelist text of a value (as a par file writes it)."
format_value(v) = format_scalar(v)
format_value(v::AbstractVector) = join(format_scalar.(v), " ")

end # module

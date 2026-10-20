# Global sea level of Spratt and Lisiecki (2016) for [bsl] sl_path: the short
# stack (sl_short, seven reconstructions, 0-430 ka, every 1 kyr on the LR04 age
# model) from FesmData SeaLevel-Spratt2016, shifted to 0 at the present (the
# stack is +8.5 m at 0 ka).
#
# Usage (from the yelmox root, fesmdata link in place):
#   julia --project=scripts/grl-paleo scripts/grl-paleo/make_sealevel_spratt2016.jl
#
# Writes input/sealevel_spratt2016.dat (time [yr, negative before present], sea level [m]),
# with one header line: FastIsostasy read_series skips the first line only (a
# later comment line leaves an unset entry).

using NCDatasets
using Printf

const SRC = "fesmdata/Records/SeaLevel/SeaLevel-Spratt2016.nc"
const OUT = "input/sealevel_spratt2016.dat"

ds = NCDataset(SRC)
t  = -Float64.(ds["time"][:])                     # years ago -> years relative to the present
sl = Float64.(coalesce.(vec(ds["sl_short"][:, :]), NaN))
close(ds)

ok = .!isnan.(sl)
t, sl = t[ok], sl[ok]
i = sortperm(t)
t, sl = t[i], sl[i]
k0 = findfirst(==(0.0), t)
k0 === nothing && error("no value at 0 ka in $SRC")
sl0 = sl[k0]
sl .-= sl0

open(OUT, "w") do io
    @printf(io, "# Time[yr]\tRSL[m] (Spratt and Lisiecki 2016 sl_short, LR04 ages, %+.2f m at 0 ka removed; scripts/grl-paleo/make_sealevel_spratt2016.jl)\n", sl0)
    for (tk, sk) in zip(t, sl)
        @printf(io, "%d\t%.4f\n", round(Int, tk), sk)
    end
end
@printf("%s: %d values, %.0f to %.0f yr, min %.1f m at %.0f yr, %.1f m at -21 kyr\n",
        OUT, length(t), t[1], t[end], minimum(sl), t[argmin(sl)], sl[findfirst(==(-21000.0), t)])

#!/usr/bin/julia

import Pkg; Pkg.activate(".")

#using IceAnalysis
import IceAnalysis

#include("/Users/robinson/models/IceAnalysis.jl/src/IceAnalysis.jl")
#using .IceAnalysis



if length(ARGS) != 1
    throw("Error: path to ensemble must be provided.")
end

path = ARGS[end];


IceAnalysis.ensemble_check(path);

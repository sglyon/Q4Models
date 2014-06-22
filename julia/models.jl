#=

Generic Model interface

@author: Spencer Lyon <sgl290@stern.nyu.edu>
@date: 06/21/2014 11:19 PM

=#

module Models

using Grid

import Base: show, getindex

## PATCH: Vectorized evaluation at multiple points for Grid
function getindex{T,R<:Real}(G::InterpIrregular{T,1}, x::AbstractVector{R})
    n = length(x)
    v = Array(T, n)
    for i = 1:n
        v[i] = getindex(G, x[i])
    end
    return v
end

#########################
#- Type: AbstractModel -#
#########################

abstract AbstractModel

function show(io::IO, mod::AbstractModel)
    # model header
    if isdefined(mod, :mod_name)
        msg = getfield(mod, :mod_name)
    else
        msg = "Generic Model"
    end
    msg *= "\n"
    msg *= "-" ^ (length(msg) - 1)

    # listing of parameters
    param_str = ""
    param_syms = names(mod)
    n = maximum(map(length, map(string, param_syms)))
    for p in param_syms
        if is(p, :Sol)
            if isdefined(mod, p)
                param_str *= "  $(rpad(string(p), n)): Computed\n"
            else
                param_str *= "  $(rpad(string(p), n)): Not yet computed\n"
            end
            continue
        end

        if is(p, :Sim)
            if isdefined(mod, p)
                f = getfield(mod, p)
                param_str *= "  $(rpad(string(p), n)): Computed ($(f.T) × $(f.N))\n"
            else
                param_str *= "  $(rpad(string(p), n)): Not yet computed\n"
            end
            continue
        end

        # don't print this twice
        is(p, :mod_name) && continue

        f = getfield(mod, p)
        if isa(f, Array)
            param_str *= "  $(rpad(string(p), n)): $(summary(f))\n"
        else
            param_str *= "  $(rpad(string(p), n)): $f\n"
        end

    end

    # put it together
    msg *= "\n" * param_str

    print(io, msg)
end


function solve!(mod::AbstractModel)
    error("This function must be implemented directly by subtypes")
end

function simulate!(mod::AbstractModel)
    error("This function must be implemented directly by subtypes")
end


############################
#- Type: AbstractSolution -#
############################

abstract AbstractSolution

##############################
#- Type: AbstractSimulation -#
##############################

abstract AbstractSimulation

####################
#- include/export -#
####################

include("ifp.jl")

export
    AbstractModel, solve!, simulate!,
    AbstractSolution,

    IFP


end  # module


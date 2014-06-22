#=

Module to hold random tools collected along the way

@author: Spencer Lyon <sgl290@stern.nyu.edu>
@date: 06/21/2014 11:39 PM

=#

module Tools

import Base: show

include("markov.jl")
include("discretize.jl")

export
    # markov
    Markov,
    simulate,

    # discretize
    tauchen,
    rouwenhorst

end  # Module

#=

Generic Model interface

@author: Spencer Lyon <sgl290@stern.nyu.edu>
@date: 06/21/2014 11:39 PM

=#

module Tools

export
    # markov
    Markov,
    simulate,

    # discretize
    tauchen,
    rouwenhorst

include("markov.jl")
include("discretize.jl")

end  # Module

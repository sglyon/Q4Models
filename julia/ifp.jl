#=

Type to define the income fluctuation problem (consumption-savings
problem)

@author : Spencer Lyon <spencer.lyon@nyu.edu>
@date : 2014-04-19 20:10:42

=#

import Tools: Markov, simulate, rouwenhorst

type IFPSol{R <: FloatingPoint} <: AbstractSolution
    c::Matrix{R}
    ap::Matrix{R}
end

type IFPSim{R <: FloatingPoint} <: AbstractSimulation
    T::Integer
    N::Integer
    burn::Integer
    c_it::Matrix{R}
    a_it::Matrix{R}
    y_it::Matrix{R}
end


type IFP{T <: FloatingPoint} <: AbstractModel
    # Model parameters
    rho::T
    beta::T
    r::T
    R::T
    gamma::T

    # Grid parameters for a
    n_a::Integer
    a_min::T
    a_max::T
    a::Vector{T}

    # Parameters for y
    n_y::Integer
    sig_eps::T
    y::Vector{T}
    P::Matrix{T}
    m::Markov

    # model type
    mod_name::String

    # Hold the solution and simulation
    Sol::IFPSol
    Sim::IFPSim

    # Incomplete initialization. Everything but Sol/Sim fields
    function IFP(rho::T, beta::T, r::T, R::T, gamma::T,
                 n_a::Integer, a_min::T, a_max::T, a::Vector{T},
                 n_y::Integer, sig_eps::T, y::Vector{T}, P::Matrix{T},
                 m::Markov)

        new(rho, beta, r, R, gamma, n_a, a_min, a_max, a,
            n_y, sig_eps, y, P, m, "Income Fluctuation Problem")
    end

    # Default constructor
    function IFP(rho::T, beta::T, r::T, R::T, gamma::T,
                 n_a::Integer, a_min::T, a_max::T, a::Vector{T},
                 n_y::Integer, sig_eps::T, y::Vector{T}, P::Matrix{T},
                 m::Markov, mod_name::String, Sol::IFPSol, Sim::IFPSim)

        new(rho, beta, r, R, gamma, n_a, a_min, a_max, a, n_y,
            sig_eps, y, P, m,  "Income Fluctuation Problem", Sol)
    end

end


function IFP{T <: FloatingPoint}(;rho::T=0.90, beta::T=0.95, r::T=0.02,
                                 gamma::T=2.0,
                                 n_a::Integer=300, a_min::T=0.0, a_max::T=100.,
                                 n_y::Integer=5, sig_eps::T=sqrt(0.06))

    # Compute a
    a::Vector{T} = linspace(a_min, a_max, n_a)

    # Compute R
    R = 1.0 + r

    # compute P, y, m
    wbar = - sig_eps^2 / (2 * (1 + rho))
    w::Vector{T}, P::Matrix{T} = rouwenhorst(n_y, rho, sig_eps, wbar)
    y = exp(w)
    m = Markov(P, y)

    return IFP{T}(rho, beta, r, R, gamma, n_a, a_min, a_max, a, n_y, sig_eps,
                  y, P, m)
end


# Define utility function and marginal utility and its inverse
# NOTE: These will all be inlined
u{T <: Real}(c::Union(T, Array{T}), γ::Float64) = c.^(1.0 - γ) ./ (1.0 - γ)
up{T <: Real}(c::Union(T, Array{T}), γ::Float64) = c.^(- γ)
up_inv{T <: Real}(b::Union(T, Array{T}), γ::Float64) = b.^(- 1.0 / γ)


#=
    Set up the initial guess for the policy function of c.

    It uses the analytical solution to the same problem with quadratic
    utility and iid shocks.
=#
function initialize_c(a::Vector{Float64}, y::Vector{Float64}, r::Float64)
    n_a, n_y = map(length, {a, y})
    c = Array(Float64, n_a, n_y)
    for i=1:n_a, j=1:n_y
        c[i, j] = r*a[i] + y[j]
    end
    return c
end


#=
    Update the guess for c using linear interpolation. This routine is
    called at the end of each loop in the main endogenous grid
    algorithm.
=#
# Updates c_til_i inplace
function update_c!(a_star::Matrix{Float64}, c_til::Matrix{Float64},
                   c_til_i::Matrix{Float64}, a_star_1::Vector{Float64},
                   a::Vector{Float64}, y::Vector{Float64}, R::Float64)

    num_bind = 0
    n_a, n_y = size(c_til)
    for j=1:n_y  # Do interpolation (step 7)
        ig = InterpIrregular(a_star[:, j], c_til[:, j], BCnearest,
                             InterpLinear)
        c_til_i[:, j] = ig[a]

        # update boundary condition on lower end
        for i=1:n_a
            if a[i] < a_star_1[j]
                num_bind += 1
                c_til_i[i, j] = R * a[i] + y[j] - a[1]
            end
        end
    end
    return num_bind
end


#=
    Apply the endogenous grid method to solve for the policy functions
    c(a, y) and a'(a, y) for a given specification of the coefficient of
    absolute risk aversion. The return value is two (n_a, n_y) matrices
    representing the optimal policy rules for each value on the grids
    for a and y.

    tol tolerance level for convergence
=#
function endog_grid!(mod::IFP; tol::Float64=1e-12, verbose=false, max_it=600)
    # Unpack parameters from model object
    n_a, n_y, a = mod.n_a, mod.n_y, mod.a

    # step 1 (discretize process for y) was done when model obj was created
    y, P, beta, gamma, r, R = (mod.y, mod.P, mod.beta, mod.gamma,
                               mod.r, mod.R)

    # Allocate memory for arrays used in iteration
    B = Array(Float64, n_a, n_y)
    c_til = Array(Float64, n_a, n_y)
    a_star = Array(Float64, n_a, n_y)
    a_star_1 = Array(Float64, n_y)
    c_til_i = Array(Float64, n_a, n_y)
    c = initialize_c(a, y, r)  # step 2

    # set up iteration parameters
    err = 1.0
    it = 0

    # Main algorithm
    while err > tol
        B[:] = beta * R .* (up(c, gamma) * P')  # step 3
        c_til[:] = up_inv(B, gamma)::Matrix{Float64}  # step 4
        a_star[:] = (c_til .+ a .- y') ./ R  # step 5
        a_star_1[:] = a_star[1, :]  # step 6
        bound = update_c!(a_star, c_til, c_til_i, a_star_1, a, y, R)  # step 7
        err = Base.maxabs(c - c_til_i)  # step 8

        # Update guess for c
        c[:] = c_til_i

        # print status, break if running too long
        it += 1
        verbose && it % 5 == 0 && println("it=$it\tbound=$bound\terr=$err")

        it > max_it && break
    end

    println("Converged. Total iterations: $it")

    # compute a'(a, y) from budget constraint using c
    ap = (R .* a) .+  y' - c

    # Attach policy rules to model object
    mod.Sol = IFPSol(c, ap)

    return nothing # don't return anything. model is updated in place
end

# fit in with AbstractModel api
solve!(mod::IFP; k...) = endog_grid!(mod; k...)


#=
    Given a set of policy functions, the Markov process for y,
    and a desired simulation length, simulate the economy for a
    specified number of periods. The return value is three Vectors: one
    for the path of consumption, one for the path of asset holdings,
    and a third for the path of the endowment process.

    a0_i is the index of the starting level of asset holdings.

    burn is the burn in period. It is set to remove 10% of simulation
    length by default.

    seed is the seed for the random number generator. The default value
    is a random number (what?!! bootstrapping the RNG?!)
=#
function simulation(mod::IFP; T::Integer=20000, a0::Number=mod.a_min,
                    N::Integer=1, burn::Integer=2000,
                    seed::Number=int(500*rand()))
    if ~isdefined(mod, :Sol)
        msg = "Sol field not found on model.\n"
        msg *= "Must call endog_grid! before simulation"
        msg *= "I'll try to do it for you"
        println(msg)
        endog_grid!(model, verbose=false)
    end

    # Pull out needed objects from model
    c::Matrix{Float64} = mod.Sol.c
    a_prime::Matrix{Float64} = mod.Sol.ap
    m::Markov = mod.m
    n_a::Integer = mod.n_a
    n_y::Integer = mod.n_y
    R::Float64 = mod.R

    # Set up data storage
    c_data = Array(Float64, T, N)
    y_data = Array(Float64, T, N)
    a_data = Array(Float64, T+1, N)

    # seed RNG to get consistent simulation of Markov process
    srand(seed)

    # Construct interpolator objects to use in computing c inside loop
    igs = [InterpIrregular(a_prime[:, j], c[:, j], BCnearest, InterpLinear)
           for j=1:n_y]

    i_t = Array(Int8, T)
    for i=1:N
        # simulate markov process
        y_t, i_t = simulate(m, T, ret_ind=true)
        y_data[:, i] = y_t

        # use initial condition to start income process
        a_data[1, i] = a0

        # simulate forward all T time periods
        for t=1:T
            @inbounds begin
            c_data[t, i] = igs[i_t[t]][a_data[t, i]]
            a_data[t+1, i] = R * a_data[t, i] + y_data[t, i] - c_data[t, i]
            end
        end
    end

    mod.Sim = IFPSim(T, N, burn, c_data[burn+1:end, :], a_data[burn+1:end, :],
               y_data[burn+1:end, :])

    return nothing
end

simulate!(mod::IFP; k...) = simulation(mod; k...)

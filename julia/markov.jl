import Base.show

type Markov
    P::Matrix{Float64}
    y::Vector{Float64}   # vector of state values
    pi_0::Union(Vector{Float64}, Vector{Int64})

    function Markov(P, y, pi_0)
        if size(P, 1) != size(P, 2)
            error("P must be a square matrix")
        end

        if size(pi_0, 1) != size(P, 1)
            error("pi_0 must have same number of elements as each row of P")
        end

        if size(y, 1) != size(P, 1)
            error("y must have same number of elements as each row of P")
        end

        new(P, y, pi_0)
    end
end


function Markov(P::Matrix{Float64}, y::Vector{Float64})
    pi_0 = zeros(size(P, 1))
    pi_0[1] = 1.0
    Markov(P, y, pi_0)
end


function Markov(P::Matrix{Float64})
    y = [1.:size(P, 1)]
    pi_0 = zeros(size(P, 1))
    pi_0[1] = 1.0
    Markov(P, y, pi_0)
end


function show(io::IO, m::Markov)
    msg = "$(length(m.y)) state Markov Chain"
    print(io, msg)
end


function simulate(mc::Markov, N::Int; ret_ind=false)
    # take cumsum along each row of transition matrix
    cs_P = cumsum(mc.P, 2)

    # Draw N random numbers on U[0, 1) and allocate memory for draw
    α = rand(N)
    draw = Array(Int16, N)

    # Compute starting value using initial distribution pi_0
    draw[1] = findfirst(α[1] .< cumsum(mc.pi_0))

    # Fill in chain using transition matrix
    for j=2:N
        draw[j] = findfirst(α[j] .< cs_P[draw[j-1], :])
    end

    if ret_ind
        return mc.y[draw], draw
    else
        return mc.y[draw]
    end
end

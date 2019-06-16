# Basic implementation of the simplex algorithm

struct Result{T}
    A        :: Matrix{T}
    b        :: Vector{T}
    c        :: Vector{T}
    feasible :: Bool
    solution :: Vector{T}
    basis    :: Vector{Int}
end

infeasible(A, b, c) = Result{Nothing}(A, b, c, false, [], [])
Result(A, b, c, vertices::Vector{T}, basis::Vector{Int}) where {T} = Result{T}(A, b, c, true, vertices, basis)


fix_simplex = false

function simplex(A::Matrix, b::Vector, c::Vector, basis::Vector{Int})
    nPlanes, nVars = size(A)
    @assert size(b) == (nPlanes,)
    @assert size(c) == (nVars,)
    entering_index = 0

    while true
        AA = A[basis, :]
        AAinv = try
            inv(AA)
        catch e
            error("Singular matrix for entering index $entering_index")
        end
        bb = b[basis]

        primal_vertex = AAinv * bb
        dual_vertex   = transpose(AAinv) * c

        i = argmin(dual_vertex)
        dual_vertex[i]
        if dual_vertex[i] >= 0
            return Result(A, b, c, primal_vertex, basis)
        end

        leaving_index = basis[i]
        direction_rhs = zeros(nPlanes)
        direction_rhs[leaving_index] = 1
        direction_rhs = direction_rhs[basis]
        edge_direction = AAinv * direction_rhs

        hyperplanes = A*edge_direction
        distances = A*primal_vertex - b

        entering_index = 0
        min_scaled_distance = Inf
        for i in 1:nPlanes
            hyperplanes[i] >= 0 && continue

            # Caution: bug!!!
            if fix_simplex
                i in basis && continue
            end

            scaled_distance = -distances[i] / hyperplanes[i]
            if scaled_distance < min_scaled_distance
                min_scaled_distance = scaled_distance
                entering_index = i
            end
        end

        if entering_index == 0
            return infeasible(A, b, c)
        end

        splice!(basis, i, [])
        push!(basis, entering_index)
    end
end

function simplex(A::Matrix, b::Vector, c::Vector)
    m, n = size(A)

    A1 = zeros(eltype(A), m+n+2, n+1)
    A1[1:m, 1:n] = A
    for i in 1:n
        A1[m+i,i] = 1
    end
    A1[1:m, n+1] = -b

    b1 = zeros(eltype(A), m+n+2)

    A1[m+n+1,n+1] = 1
    A1[m+n+2,n+1] = -1

    b1[m+n+2, 1] = -1

    c1 = zeros(n+1)
    c1[n+1] = -1

    basis1 = collect(m+1:m+n+1)
    res = simplex(A1, b1, c1, basis1)

    basis2 = res.basis

    i = findfirst(isequal(m+n+2), basis2)
    i === nothing && return infeasible(A, b, c)
    splice!(basis2, i, [])
    basis2

    A2 = A1[1:m+n,1:n]
    b2 = -A1[1:m+n, n+1]
    simplex(A2, b2, c, basis2)
end



# Helper functions

load(T, fname) = open(fname, "r") do f
    siz = readline(f) |> Meta.parse
    A = T(undef, siz.args...)
    for i in eachindex(A)
        A[i] = parse(Float64, readline(f))
    end
    A
end

using StochasticArithmetic
StochasticArithmetic.value(x::Rational) = Float64(x)

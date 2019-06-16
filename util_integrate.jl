using Plots
using BenchmarkTools

function test(integrate_err, T)
    err = []
    nmsh = Int[]

    n = 10
    while n < 10_000
        push!(nmsh, n)
        push!(err,  integrate_err(n, T))

        n = Int(round(n*1.2))
    end

    p = plot(nmsh, [e.ref for e in err],
         label="ref",
         xaxis=("N",   :log10),
         yaxis=("err", :log10),
         ylims=(eps(T),1))

    for col in fieldnames(typeof(err[1]))[2:end]
        scatter!(nmsh, [getproperty(e, col) for e in err],
             label=col)
    end

    p
end







"""
Rational number, with numerator and denominator of type T
"""
struct MyRational{T<:Integer} <: Real
    _num   :: T
    _denom :: T

    function MyRational{T}(n,d) where {T<:Integer}
        d==0 && error("zero denominator")
        new(n,d)
    end
end

MyRational(num::T, denom::T) where {T} = MyRational{T}(num, denom)
MyRational{T}(i::Int) where {T} = MyRational{T}(i, 1)


num(a::MyRational)   = a._num
denom(a::MyRational) = a._denom

function my_gcd(a, b)
    r = 1
    while r != 0
        r = mod(a, b)
        (a, b) = (b, r)
    end
    a
end

function simplify_rational(num, denom)
    if denom < 0
        (num, denom) = (-num, -denom)
    end
    @assert denom > 0
    d = my_gcd(num, denom)
    MyRational(div(num, d), div(denom, d))
end

Base.promote_rule(::Type{MyRational{T}}, ::Type{Int}) where {T} = MyRational{T}

import Base.*
function *(a::MyRational, b::MyRational)
    simplify_rational(num(a)*num(b), denom(a)*denom(b))
end

import Base./
function /(a::MyRational, b::MyRational)
    simplify_rational(num(a)*denom(b), denom(a)*num(b))
end

import Base.+
function +(a::MyRational, b::MyRational)
    simplify_rational(num(a)*denom(b)+num(b)*denom(a), denom(a)*denom(b))
end

Base.:-(a::MyRational) = MyRational(-num(a), denom(a))
Base.:-(a::MyRational, b::MyRational) = a + (-b)

import Base.<
function <(a::MyRational, b::MyRational)
    num(b) == 0 && return num(a) < 0
    r = a / b
    num(r) < denom(r)
end

import Base.==
function ==(a::MyRational, b::MyRational)
    num(a) == num(b) && denom(a) == denom(b)
end

Base.show(io::IO, a::MyRational{T}) where {T} = print(io, num(a), "/", denom(a))

to_float(x::MyRational) = num(x)/denom(x)

nothing

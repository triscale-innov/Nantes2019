function read_data(fname)
    open(fname) do f
        n, cond, ref = split(readline(f), " ")
        cond, ref = parse.(Float64, (cond, ref))
        x = parse.(Float64, eachline(f))
        cond, ref, x
    end
end

const TESTDIR = "dataSum"

function test(fun, typ)
    cond = Float64[]
    nsig = Float64[]
    map(readdir(TESTDIR)) do fname
        χ, ρ, x = read_data(joinpath(TESTDIR, fname))
        σ = reliable_digits(x, fun, typ(ρ))

        push!(cond, χ)
        push!(nsig, σ)
    end

    label = string(nameof(fun)) * ", " * string(typ)
    scatter!(cond, nsig, label=label)
end

using Plots

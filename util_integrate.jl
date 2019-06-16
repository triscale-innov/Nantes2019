using Plots

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

nothing

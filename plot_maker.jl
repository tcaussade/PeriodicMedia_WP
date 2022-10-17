using Plots
using DelimitedFiles

x = Vector{Float64}(undef,0)
open("A_windowconvergence_ppw1_dim1_k8.8.txt") do io
    while !eof(io)
        xt    = readline(io)
        append!(x,parse(Float64,xt))
    end
end

y = Vector{Float64}(undef,0)
open("eb_windowconvergence_ppw1_dim1_k8.8.txt") do io
    while !eof(io)
        yt    = readline(io)
        append!(y,parse(Float64,yt))
    end
end

p1 = plot(x,log10.(y), xlabel = "A/λ", ylabel = "Error", title = "semilog")
p2 = plot(log10.(x), log10.(y), xlabel = "A/λ", ylabel = "Error", title = "loglog")
plot(p1,p2, legend = false)
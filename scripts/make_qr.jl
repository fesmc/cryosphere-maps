# One-off: write the QR code of the website URL to data/qr_site.txt (rows of 0/1,
# 1 = dark module) for the posters. QRCoders is installed in a temporary
# environment because it holds back CairoMakie's dependencies.
# Usage: julia scripts/make_qr.jl [url]
using Pkg
Pkg.activate(; temp=true); Pkg.add("QRCoders"; io=devnull)
using QRCoders

url = isempty(ARGS) ? "https://fesmc.github.io/cryosphere-maps/" : ARGS[1]
m = qrcode(url; eclevel=Medium())
out = joinpath(@__DIR__, "..", "data", "qr_site.txt")
open(out, "w") do io
    println(io, "# ", url)
    for r in eachrow(m)
        println(io, join(Int.(r)))
    end
end
println("wrote ", out, " (", size(m, 1), "×", size(m, 2), ")")

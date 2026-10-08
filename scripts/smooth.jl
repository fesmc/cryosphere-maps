# Gaussian smoothing of grids (shared by the plot and prepare scripts).

"Separable Gaussian smoothing (sigma in grid cells), NaN-aware."
function gauss_smooth(z, σ)
    r = ceil(Int, 3σ)
    k = [exp(-0.5*(i/σ)^2) for i in -r:r]
    function pass(a, dim)
        b = similar(a); nx, ny = size(a)
        for j in 1:ny, i in 1:nx
            s = 0.0; w = 0.0
            for (q, kq) in zip(-r:r, k)
                ii, jj = dim == 1 ? (i+q, j) : (i, j+q)
                (1 <= ii <= nx && 1 <= jj <= ny) || continue
                v = a[ii, jj]; isnan(v) && continue
                s += kq*v; w += kq
            end
            b[i, j] = w > 0 ? s/w : NaN
        end
        return b
    end
    return pass(pass(z, 1), 2)
end

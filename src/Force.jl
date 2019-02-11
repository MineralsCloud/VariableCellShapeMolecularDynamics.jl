"""
# module Force



# Examples

```jldoctest
julia>
```
"""
module Force

using Distances: euclidean

using VCSMD.Geometry

export lennard_jones_scalar,
    lennard_jones_vector

function lennard_jones_scalar(ε::Float64, σ::Float64)::Function
    function (r::Float64)
        x = σ / r
        24 * ε / σ * (2 * x^13 - x^7)
    end
end

function lennard_jones_vector(ε::Float64, σ::Float64)::Function
    f = lennard_jones_scalar(ε, σ)

    function (ri::Point{N}, rj::Point{N}) where {N}
        r = euclidean(rj, ri)
        (rj - ri) / r * f(r)
    end
end

end
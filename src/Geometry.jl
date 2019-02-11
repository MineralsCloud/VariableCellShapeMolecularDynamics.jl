"""
# module Geometry



# Examples

```jldoctest
julia>
```
"""
module Geometry

using LinearAlgebra: det

using StaticArrays: FieldVector, SMatrix

export Point,
    Point2D, Point3D,
    Cell,
    volume,
    metric_tensor,
    inverse_metric_tensor,
    metric_tensor_velocity,
    reciprocalcell

abstract type Point{N, Float64} <: FieldVector{N, Float64} end

struct Point2D <: Point{2, Float64}
    x::Float64
    y::Float64
end

struct Point3D <: Point{3, Float64}
    x::Float64
    y::Float64
    z::Float64
end

const Cell = SMatrix{3, 3, Float64}

volume(c::Cell)::Float64 = (abs ∘ det)(c)

metric_tensor(c::Cell) = transpose(c) * c

inverse_metric_tensor(c::Cell) = (inv ∘ metric_tensor)(c)

metric_tensor_velocity(c::Cell, c_derivative::Cell) = transpose(c_derivative) * c + c * transpose(c_derivative)

reciprocalcell(c::Cell) = 2π * (transpose ∘ inv)(c)

end

"""
# module Geometry



# Examples

```jldoctest
julia>
```
"""
module Geometry

using StaticArrays: FieldVector

export Point,
    Point2D, Point3D

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

end
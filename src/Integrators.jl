"""
# module Integrators



# Examples

```jldoctest
julia>
```
"""
module Integrators

export velocity_verlet

function velocity_verlet(x::Float64, v::Float64, f::Function, dt::Float64, atomic_mass::Float64)
    nextx = x + v * dt + f(x) / (2 * atomic_mass) * dt^2
    nextv = v + (f(nextx) + f(x)) / (2 * atomic_mass) * dt
    nextx, nextv
end

function velocity_verlet(x::Float64, v::Float64, f::Function, dt::Float64, atomic_mass::Float64, steps_amount::Int)
    for i in 1:steps_amount
        x, v = velocity_verlet(x, v, f, dt, atomic_mass)
    end
    x, v
end

end

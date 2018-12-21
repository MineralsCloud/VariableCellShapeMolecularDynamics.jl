"""
# module Thermostat



# Examples

```jldoctest
julia>
```
"""
module Thermostat

using Distributions: Normal
using StaticArrays: SVector, SMatrix, MVector, MMatrix

export andersen_heat_bath,
    homogeneous_atomic_gas_energy,
    perform_stochastic_collision,
    velocity_after_collision

const Boltzmann = 1

"""
    andersen_heat_bath(energy, atoms_amount[, dimension])



# Arguments
- `energy::Float64`: the total energy of atoms of an element.
- `atoms_amount::Int`: the number of atoms of an element.
- `dimension::Int=3`: the dimensionality of the problem.
"""
function andersen_heat_bath(energy::Float64, atoms_amount::Int, dimension::Int = 3)::Normal{Float64}
    instantaneous_temperature = energy / (Boltzmann * dimension * atoms_amount)
    σ = sqrt(instantaneous_temperature)  # Standard deviation of temperature
    Normal(0, σ)
end

"""
    homogeneous_atomic_gas_energy(atomic_mass, velocities)

# e.g., v_i = v_{ix}^2 + v_{iy}^2 + v_{iz}^2 if `D` is `3`
# m_i v_i^2 for the atoms of an element

# Arguments
- `atomic_mass::Float64`:
- `velocities::SMatrix{N, D, Float64}`: `N` is the number of atoms of an element, `D` is the the
  dimensionality of the problem.
"""
function homogeneous_atomic_gas_energy(atomic_mass::Float64, v::SMatrix{N, D, Float64})::Float64 where {N, D}
    atomic_mass * sum(v .^ 2)
end

"""
    perform_stochastic_collision(energy, ν, dt[, dimension])

Perform one stochastic collision with the heat bath.

# Arguments
- `energy::Float64`:
- `ν::Float64`: `N` is the number of atoms of an element, `D` is the the
  dimensionality of the problem.
- `dt::Float64`:
- `dimension::Int=3`: the dimensionality of the problem.
"""
function perform_stochastic_collision(energy::Float64, ν::Float64, dt::Float64, dimension::Int = 3)::Vector{Tuple{Int, Float64}} where {N}
    normal_distribution_random_number_generator = andersen_heat_bath(energy, N, dimension)
    probability = ν * dt  # The probability that a particle is selected in a time step of length Δt
    [rand() < probability && (i, rand(normal_distribution_random_number_generator)) for i in 1:N]
end

"""
"""
function velocity_after_collision(atomic_mass::Float64, v::SMatrix{N, D, Float64}, ν::Float64, dt::Float64)::SMatrix{N, D, Float64} where {N, D}
    energy = homogeneous_atomic_gas_energy(atomic_mass, v)
    v = MMatrix(v)
    for (i, velocity) in perform_stochastic_collision(energy, ν, dt, D)
        v[i, :] .= velocity / sqrt(D)
    end
    v
end

end

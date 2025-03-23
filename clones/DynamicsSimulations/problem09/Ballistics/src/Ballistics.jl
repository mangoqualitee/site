module Ballistics

using DifferentialEquations
using UnPack
using LinearAlgebra
using GLMakie

include("Parameters.jl")
include("Physics.jl")
include("Visualization.jl")
include("Benchmarks.jl")

# problem setup
x0, y0 = 0.0, 0.0
r0 = [x0; y0]
speed0 = 1
launchangle = pi / 4
vx0, vy0 = speed0 * [cos(launchangle); sin(launchangle)]
v0 = [vx0; vy0]
u0 = [r0; v0]

t0 = 0.0
tend = 1.0
tspan = (t0, tend)
p = Parameters.Param(m = 1, g = 1, c = 1)
prob = ODEProblem(Physics.ballistic!, u0, tspan, p)

function analytical_solve(t)
    m, g, c = p.mass, p.gravity, p.viscosity

    x =
        x0 +
        vx0 * ((-m / c) * exp((c / m) * (t0))) * (exp((-c / m) * t) - exp((-c / m) * t0))
    y =
        y0 +
        (
            (vy0 + (m * g / c)) *
            (-m / c) *
            (exp((c / m) * t0)) *
            (exp((-c / m) * t) - exp((-c / m) * t0))
        ) +
        ((-m * g / c) * (t - t0))
    vx = (vx0) * (exp((-c / m) * (t - t0)))
    vy = ((vy0) + (m * g / c)) * exp((-c / m) * (t - t0)) - (m * g / c)

    u = zeros((4,))
    u[1] = x
    u[2] = y
    u[3] = vx
    u[4] = vy

    return u
end

# solver
# Δh = 0.1
# sol_numeric = solve(prob, Midpoint(), dt=Δh, adaptive=false)
# sol_analytical = analytical_solve.(t0:Δh:tend)

# visualize
# Visualization.plot_trajectory(sol_numeric.u)
# Visualization.plot_both_trajectory(sol_numeric.u, sol_analytical)
# Visualization.plot_both_position(sol_numeric.u, sol_analytical)
# Visualization.plot_both_velocity(sol_numeric.u, sol_analytical)


# log_steps, log_errors_dict = Benchmarks.abs_error_vs_step_sizes(prob, analytical_solve)
# Visualization.plot_abs_steps_vs_error(log_steps, log_errors_dict)

# time_histories, error_histories = Benchmarks.error_vs_time(prob, analytical_solve)
# Visualization.plot_abs_steps_vs_time(time_histories, error_histories)

# (xs, ys, zs) = Benchmarks.tolerances_vs_error(prob, analytical_solve)
# Visualization.plot_tolerances_vs_error(xs, ys, zs)

end # module Ballistics

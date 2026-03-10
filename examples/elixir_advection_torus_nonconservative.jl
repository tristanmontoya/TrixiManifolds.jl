###############################################################################
# Torus advection test (covariant nonconservative form)
#
# This elixir constructs a periodic torus mesh, computes torus covariant metric
# terms, runs covariant linear advection in nonconservative form with a Gaussian
# scalar bump and a tangential velocity field, and converts output to VTK.
# Note: this example is unstable and is intended for demonstration purposes only.
###############################################################################
using OrdinaryDiffEq
using Trixi, TrixiAtmo, Trixi2Vtk
using TrixiManifolds

# Torus geometry and mesh resolution
const major_radius = 1.0
const minor_radius = 0.5 * major_radius
const cells_per_major_angle = 32
const cells_per_minor_angle = 12
const output_directory = normpath(@__DIR__, "output", "advection_torus_nonconservative")
const rotation_period = 1.0
const num_revolutions = 4
const final_time = num_revolutions * rotation_period
const save_interval = 0.02

# Gaussian scalar bump and solid-body rotation around the torus with period T = 1
const theta_angular_speed = 2pi / rotation_period
const theta_bump = pi / 2
const phi_bump = 0.0
const gaussian_height = 1.0
const gaussian_width = 5.0

# Initial condition for scalar Gaussian bump in global Cartesian coordinates
@inline function initial_condition_transport(x, t, equations)
    bump_centre = torus2cartesian(theta_bump, phi_bump, major_radius, minor_radius)
    x1, x2, x3 = x
    dx1 = x1 - bump_centre[1]
    dx2 = x2 - bump_centre[2]
    dx3 = x3 - bump_centre[3]
    h = gaussian_height * exp(-gaussian_width * (dx1^2 + dx2^2 + dx3^2) / major_radius^2)
    return SVector(h)
end

# Global Cartesian solid-body velocity on the torus
@inline function velocity_global_cartesian(x)
    x1, x2, _x3 = x
    return theta_angular_speed * SVector(-x2, x1, 0.0)
end

# Scalar advection coefficients in local contravariant coordinates
@inline function advection_coefficient_matrix(x, aux_vars, orientation::Integer, equations)
    velocity_contravariant = TrixiAtmo.basis_contravariant(aux_vars, equations) *
                             velocity_global_cartesian(x)
    return SMatrix{1, 1}(velocity_contravariant[orientation])
end

# Conservative coefficients are zero in this nonconservative example
@inline function conservative_coefficient_matrix(x, aux_vars, orientation::Integer,
                                                 equations)
    return SMatrix{1, 1}(zero(eltype(aux_vars)))
end

# Nonconservative coefficients carry the advection operator
@inline function nonconservative_coefficient_matrix(x, aux_vars, orientation::Integer,
                                                    equations)
    return advection_coefficient_matrix(x, aux_vars, orientation, equations)
end

# Clear stale output files from previous runs to avoid mixing datasets.
if isdir(output_directory)
    rm(output_directory; recursive = true, force = true)
end
mkpath(output_directory)

# Covariant 1x1 linear system in nonconservative form with B^j ∂_j u
equations = CovariantLinearSystem2D(1, conservative_coefficient_matrix,
                                    nonconservative_coefficient_matrix,
                                    global_coordinate_system = GlobalCartesianCoordinates())

# DGSEM solver with nonconservative surface and volume terms
volume_flux = (flux_central, flux_nonconservative)
surface_flux = (flux_lax_friedrichs, flux_nonconservative)
solver = DGSEM(polydeg = 3, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

# Create periodic torus mesh in 3D
mesh = P4estMeshTorus2D(cells_per_major_angle, cells_per_minor_angle,
                        major_radius, minor_radius,
                        polydeg = Trixi.polydeg(solver))

# Transform initial condition to the internal conservative variables
initial_condition_transformed = transform_initial_condition(initial_condition_transport,
                                                            equations)

# Set up semidiscretization and run window
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_transformed, solver,
                                    boundary_conditions = boundary_condition_periodic,
                                    metric_terms = MetricTermsCovariantTorus(major_radius,
                                                                             minor_radius,
                                                                             cells_per_major_angle,
                                                                             cells_per_minor_angle))
ode = semidiscretize(semi, (0.0, final_time))
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi, interval = 10, save_analysis = true,
                                     extra_analysis_errors = (:conservation_error,))
save_solution = SaveSolutionCallback(dt = save_interval,
                                     solution_variables = contravariant2global,
                                     save_initial_solution = true,
                                     output_directory = output_directory)
stepsize_callback = StepsizeCallback(cfl = 0.7)
callbacks = CallbackSet(summary_callback, analysis_callback, save_solution,
                        stepsize_callback)

# Solve the ODE using OrdinaryDiffEq.jl with low-storage Runge-Kutta method
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0e-3, save_everystep = false, callback = callbacks)

# Convert to VTU for visualization
trixi2vtk(output_directory * "/solution_*.h5", output_directory = output_directory)

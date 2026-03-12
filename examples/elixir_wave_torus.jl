###############################################################################
# Torus wave test (covariant mixed conservative/nonconservative form)
#
# This elixir constructs a periodic torus mesh, computes torus covariant metric
# terms, and solves a generic linear wave system with state
# q = (p, u¹, u²), where:
# - p is the scalar wave variable
# - u¹ and u² are contravariant velocity components
#
# The equations are implemented as
#   ∂ₜq + (1/J)∂ⱼ(JAʲq) + Bʲ∂ⱼq = 0,
# using conservative Aʲ terms for div(u) and nonconservative Bʲ terms for
# c² grad(p)
###############################################################################
using OrdinaryDiffEq
using Trixi, TrixiAtmo, Trixi2Vtk
using TrixiManifolds

# Rognes-like torus geometry and mesh resolution
const major_radius = 1.0
const minor_radius = 0.25
const cells_per_major_angle = 32
const cells_per_minor_angle = 12
const output_directory = normpath(@__DIR__, "output", "wave_torus")
const final_time = 10.0
const save_interval = 0.02

# Generic wave-speed parameter c and c² (Rognes benchmark uses c² = 1)
const wave_speed_squared = 1.0

# Initial condition evaluated at Cartesian node x = (x₁, x₂, x₃):
# q = (p, u¹, u²), with torus angular components mapped to contravariant components
@inline function initial_condition_wave(x, t, aux_vars, equations)
    x1, x2, x3 = x
    ring_radius = sqrt(x1^2 + x2^2)
    theta = atan(x2, x1)
    phi = atan(x3, ring_radius - major_radius)

    shape = (1 - 0.25 * cos(phi)) * cos(phi)
    p = 0.1 * shape
    u_theta = 0.0
    u_phi = -0.4 * shape

    sin_theta, cos_theta = sincos(theta)
    sin_phi, cos_phi = sincos(phi)
    basis_theta = SVector(-ring_radius * sin_theta, ring_radius * cos_theta, 0.0)
    basis_phi = SVector(-minor_radius * sin_phi * cos_theta,
                        -minor_radius * sin_phi * sin_theta,
                        minor_radius * cos_phi)

    velocity_global = u_theta * basis_theta + u_phi * basis_phi
    velocity_contravariant = TrixiAtmo.basis_contravariant(aux_vars, equations) *
                             velocity_global

    return SVector(p, velocity_contravariant[1], velocity_contravariant[2])
end

# Conservative coefficient matrices A¹ and A²
# A¹[1,2] = 1, A²[1,3] = 1; all other entries are zero
@inline function conservative_coefficient_matrix(x, aux_vars, orientation::Integer,
                                                 equations)
    zero_aux = zero(eltype(aux_vars))
    one_aux = one(eltype(aux_vars))
    if orientation == 1
        return SMatrix{3, 3}(zero_aux, zero_aux, zero_aux,
                             one_aux, zero_aux, zero_aux,
                             zero_aux, zero_aux, zero_aux)
    else
        return SMatrix{3, 3}(zero_aux, zero_aux, zero_aux,
                             zero_aux, zero_aux, zero_aux,
                             one_aux, zero_aux, zero_aux)
    end
end

# Nonconservative matrices Bʲ for c² grad(p) in momentum equations:
# Bʲ[2,1] = c²G¹ʲ, Bʲ[3,1] = c²G²ʲ
@inline function nonconservative_coefficient_matrix(x, aux_vars, orientation::Integer,
                                                    equations)
    zero_aux = zero(eltype(aux_vars))
    c2 = convert(eltype(aux_vars), wave_speed_squared)
    Gcon = TrixiAtmo.metric_contravariant(aux_vars, equations)

    if orientation == 1
        b21 = c2 * Gcon[1, 1]
        b31 = c2 * Gcon[2, 1]
    else
        b21 = c2 * Gcon[1, 2]
        b31 = c2 * Gcon[2, 2]
    end

    return SMatrix{3, 3}(zero_aux, b21, b31,
                         zero_aux, zero_aux, zero_aux,
                         zero_aux, zero_aux, zero_aux)
end

# Clear stale output files from previous runs to avoid mixing datasets.
if isdir(output_directory)
    rm(output_directory; recursive = true, force = true)
end
mkpath(output_directory)

# Covariant 3×3 linear wave system in mixed conservative/nonconservative form
equations = CovariantLinearSystem2D(3, conservative_coefficient_matrix,
                                    nonconservative_coefficient_matrix,
                                    global_coordinate_system = GlobalCartesianCoordinates())

# DGSEM solver with conservative and nonconservative volume/surface contributions
volume_flux = (flux_central, flux_nonconservative)
surface_flux = (flux_lax_friedrichs, flux_nonconservative)
solver = DGSEM(polydeg = 3, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

# Create periodic torus mesh in 3D
mesh = P4estMeshTorus2D(cells_per_major_angle, cells_per_minor_angle,
                        major_radius, minor_radius,
                        polydeg = Trixi.polydeg(solver))

# Set up semidiscretization and run window
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_wave, solver,
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

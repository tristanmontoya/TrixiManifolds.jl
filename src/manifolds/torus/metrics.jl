@doc raw"""
    MetricTermsCovariantTorus(major_radius, minor_radius,
                              trees_per_major_angle, trees_per_minor_angle)

Metric-term descriptor for a torus surface embedded in Cartesian 3D coordinates.
Given torus radii `R = major_radius` and `r = minor_radius`, the manifold is
parameterized by
```math
\begin{aligned}
x &= (R + r\cos\phi)\cos\theta,\\
y &= (R + r\cos\phi)\sin\theta,\\
z &= r\sin\phi.
\end{aligned}
```
The integer counts ``trees_per_major_angle = N_\theta`` and
``trees_per_minor_angle = N_\phi`` set the logical angle spacing used by the
element-to-torus map. They are required to compute basis vectors with respect to
local reference coordinates ``(\xi_1,\xi_2)`` using the correct scale factors
``d\theta/d\xi_1 = \pi/N_\theta`` and ``d\phi/d\xi_2 = \pi/N_\phi``.
"""
struct MetricTermsCovariantTorus{RealT <: Real}
    major_radius::RealT
    minor_radius::RealT
    trees_per_major_angle::Int
    trees_per_minor_angle::Int
end

function MetricTermsCovariantTorus(major_radius::Real,
                                   minor_radius::Real)
    return MetricTermsCovariantTorus(major_radius, minor_radius, 1, 1)
end

function MetricTermsCovariantTorus(major_radius::Real,
                                   minor_radius::Real,
                                   trees_per_major_angle::Integer,
                                   trees_per_minor_angle::Integer)
    RealT = promote_type(typeof(major_radius), typeof(minor_radius))
    @assert trees_per_major_angle > 0
    @assert trees_per_minor_angle > 0
    return MetricTermsCovariantTorus{RealT}(RealT(major_radius), RealT(minor_radius),
                                            Int(trees_per_major_angle),
                                            Int(trees_per_minor_angle))
end

# Initialize auxiliary variables for covariant equations on a torus in Cartesian 3D.
function TrixiAtmo.init_auxiliary_node_variables!(auxiliary_variables,
                                                  mesh::P4estMesh{2, 3},
                                                  equations::TrixiAtmo.AbstractCovariantEquations{2,
                                                                                                  3},
                                                  dg,
                                                  elements,
                                                  metric_terms::MetricTermsCovariantTorus,
                                                  bottom_topography)
    @assert equations.global_coordinate_system isa TrixiAtmo.GlobalCartesianCoordinates

    (; node_coordinates) = elements
    (; aux_node_vars) = auxiliary_variables

    major_radius = metric_terms.major_radius
    minor_radius = metric_terms.minor_radius
    trees_per_major_angle = metric_terms.trees_per_major_angle
    trees_per_minor_angle = metric_terms.trees_per_minor_angle

    @assert major_radius > zero(major_radius)
    @assert minor_radius > zero(minor_radius)
    @assert trees_per_major_angle > 0
    @assert trees_per_minor_angle > 0

    n_aux = TrixiAtmo.n_aux_node_vars(equations)
    @assert n_aux >= 26

    Trixi.@threaded for element in 1:Trixi.ncells(mesh)
        for j in eachnode(dg), i in eachnode(dg)
            x_node = Trixi.get_node_coords(node_coordinates, equations, dg, i, j, element)
            x1, x2, x3 = x_node

            # Recover torus angles from Cartesian coordinates using the same convention
            # as the (ξ₁, ξ₂) -> (θ, φ) -> (x, y, z) mapping in mesh.jl:
            # θ = 0 on +x, increasing toward +y around +z
            # φ = 0 at the outer equator (away from z-axis)
            theta = atan(x2, x1)
            rho = sqrt(x1^2 + x2^2)
            phi = atan(x3, rho - major_radius)

            s_theta, c_theta = sincos(theta)
            s_phi, c_phi = sincos(phi)

            ring_radius = major_radius + minor_radius * c_phi
            dtheta_dxi1 = convert(eltype(x_node), pi / trees_per_major_angle)
            dphi_dxi2 = convert(eltype(x_node), pi / trees_per_minor_angle)

            # Covariant basis vectors for x(θ, φ) with
            # dθ/dξ₁ = π/Nθ and dφ/dξ₂ = π/Nφ.
            dxdtheta = SVector(-ring_radius * s_theta,
                               ring_radius * c_theta,
                               zero(eltype(x_node)))
            dxdphi = SVector(-minor_radius * s_phi * c_theta,
                             -minor_radius * s_phi * s_theta,
                             minor_radius * c_phi)
            basis_covariant = SMatrix{3, 2}(dtheta_dxi1 * dxdtheta[1],
                                            dtheta_dxi1 * dxdtheta[2],
                                            dtheta_dxi1 * dxdtheta[3],
                                            dphi_dxi2 * dxdphi[1],
                                            dphi_dxi2 * dxdphi[2],
                                            dphi_dxi2 * dxdphi[3])
            aux_node_vars[1:6, i, j, element] = SVector(basis_covariant)

            metric_covariant = basis_covariant' * basis_covariant
            metric_contravariant = inv(metric_covariant)
            basis_contravariant = metric_contravariant * basis_covariant'

            aux_node_vars[7:12, i, j, element] = SVector(basis_contravariant)
            aux_node_vars[13, i, j, element] = sqrt(det(metric_covariant))
            aux_node_vars[14:16, i, j, element] = SVector(metric_covariant[1, 1],
                                                          metric_covariant[1, 2],
                                                          metric_covariant[2, 2])
            aux_node_vars[17:19, i, j, element] = SVector(metric_contravariant[1, 1],
                                                          metric_contravariant[1, 2],
                                                          metric_contravariant[2, 2])

            if !isnothing(bottom_topography)
                aux_node_vars[20, i, j, element] = bottom_topography(x_node)
            else
                aux_node_vars[20, i, j, element] = zero(eltype(aux_node_vars))
            end

            # Exact derivatives of metric components in (ξ₁, ξ₂) coordinates.
            dg11_dxi1 = zero(eltype(aux_node_vars))
            dg22_dxi1 = zero(eltype(aux_node_vars))
            dg22_dxi2 = zero(eltype(aux_node_vars))
            dg12_dxi1 = zero(eltype(aux_node_vars))
            dg12_dxi2 = zero(eltype(aux_node_vars))
            dg11_dxi2 = -2 * dtheta_dxi1^2 * minor_radius * ring_radius * s_phi * dphi_dxi2

            gcon = SMatrix{2, 2}(aux_node_vars[17, i, j, element],
                                 aux_node_vars[18, i, j, element],
                                 aux_node_vars[18, i, j, element],
                                 aux_node_vars[19, i, j, element])
            dGdxi1 = SVector(dg11_dxi1, dg12_dxi1, dg22_dxi1)
            dGdxi2 = SVector(dg11_dxi2, dg12_dxi2, dg22_dxi2)
            aux_node_vars[21:26, i, j, element] = calc_christoffel_symbols_torus(dGdxi1,
                                                                                 dGdxi2,
                                                                                 gcon)

            if n_aux > 26
                aux_node_vars[27:n_aux, i, j, element] .= zero(eltype(aux_node_vars))
            end
        end
    end

    return nothing
end

# Christoffel symbols of the second kind Γᵢⱼᵏ for 2D metric data packed as
# (dG11, dG12, dG22) derivatives and contravariant metric tensor G^{ij}.
@inline function calc_christoffel_symbols_torus(dGdxi1, dGdxi2, Gcon)
    dG11dxi1, dG12dxi1, dG22dxi1 = dGdxi1
    dG11dxi2, dG12dxi2, dG22dxi2 = dGdxi2

    gamma_1 = SMatrix{2, 2}(0.5 * dG11dxi1,
                            0.5 * dG11dxi2,
                            0.5 * dG11dxi2,
                            dG12dxi2 - 0.5 * dG22dxi1)
    gamma_2 = SMatrix{2, 2}(dG12dxi1 - 0.5 * dG11dxi2,
                            0.5 * dG22dxi1,
                            0.5 * dG22dxi1,
                            0.5 * dG22dxi2)

    return SVector(Gcon[1, 1] * gamma_1[1, 1] + Gcon[1, 2] * gamma_2[1, 1],
                   Gcon[1, 1] * gamma_1[1, 2] + Gcon[1, 2] * gamma_2[1, 2],
                   Gcon[1, 1] * gamma_1[2, 2] + Gcon[1, 2] * gamma_2[2, 2],
                   Gcon[2, 1] * gamma_1[1, 1] + Gcon[2, 2] * gamma_2[1, 1],
                   Gcon[2, 1] * gamma_1[1, 2] + Gcon[2, 2] * gamma_2[1, 2],
                   Gcon[2, 1] * gamma_1[2, 2] + Gcon[2, 2] * gamma_2[2, 2])
end

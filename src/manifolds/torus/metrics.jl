@doc raw"""
    MetricTermsCovariantTorus(major_radius, minor_radius,
                              trees_per_major_angle, trees_per_minor_angle)

Metric-term descriptor for a torus surface embedded in Cartesian 3D coordinates.
The torus map is
```math
\begin{aligned}
x &= (R + r\cos\phi)\cos\theta,\\
y &= (R + r\cos\phi)\sin\theta,\\
z &= r\sin\phi,
\end{aligned}
```
with `R = major_radius`, `r = minor_radius`.

`trees_per_major_angle = N_theta` and `trees_per_minor_angle = N_phi` are required
to map p4est tree indices to global torus angles.
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

@inline function torus_tree_indices(tree_id::Integer, trees_per_major_angle::Integer)
    # Tree ids are stored in linearized (k_θ, k_φ) order
    k_theta = Int(mod(tree_id - 1, trees_per_major_angle) + 1)
    k_phi = Int((tree_id - 1) ÷ trees_per_major_angle + 1)
    return k_theta, k_phi
end

function calc_torus_element_map_parameters(mesh::P4estMesh{2, 3},
                                           metric_terms::MetricTermsCovariantTorus,
                                           ::Type{RealT}) where {RealT <: Real}
    nelements = Trixi.ncells(mesh)
    theta_origin = Vector{RealT}(undef, nelements)
    phi_origin = Vector{RealT}(undef, nelements)
    dtheta_dxi1 = Vector{RealT}(undef, nelements)
    dphi_dxi2 = Vector{RealT}(undef, nelements)

    n_trees_expected = metric_terms.trees_per_major_angle * metric_terms.trees_per_minor_angle
    @assert Trixi.ntrees(mesh) == n_trees_expected

    two_pi = convert(RealT, 2pi)
    pi_rt = convert(RealT, pi)
    inv_p4est_root_len = inv(convert(RealT, 1 << Trixi.P4EST_MAXLEVEL))

    trees = Trixi.unsafe_wrap_sc(Trixi.p4est_tree_t, mesh.p4est.trees)
    for tree_id in eachindex(trees)
        k_theta, k_phi = torus_tree_indices(tree_id, metric_terms.trees_per_major_angle)
        @assert k_phi <= metric_terms.trees_per_minor_angle

        tree = trees[tree_id]
        quadrants = Trixi.unsafe_wrap_sc(Trixi.p4est_quadrant_t, tree.quadrants)
        tree_offset = Int(tree.quadrants_offset)

        for local_quadrant_id in eachindex(quadrants)
            quad = quadrants[local_quadrant_id]
            element = tree_offset + local_quadrant_id

            scale = p4est_quadrant_reference_scale(quad.level, RealT)

            # Global torus angles on this element are affine in local reference coordinates:
            # θ = θ_origin + (dθ/dξ₁) ⋅ (ξ₁ + 1)
            # φ = φ_origin + (dφ/dξ₂) ⋅ (ξ₂ + 1)
            # `quad.x` and `quad.y` are p4est integer coordinates in [0, 2 to the power `P4EST_MAXLEVEL`)
            theta_origin[element] = two_pi *
                                    ((k_theta - 1) + quad.x * inv_p4est_root_len) /
                                    metric_terms.trees_per_major_angle
            phi_origin[element] = two_pi *
                                  ((k_phi - 1) + quad.y * inv_p4est_root_len) /
                                  metric_terms.trees_per_minor_angle
            dtheta_dxi1[element] = pi_rt * scale / metric_terms.trees_per_major_angle
            dphi_dxi2[element] = pi_rt * scale / metric_terms.trees_per_minor_angle
        end
    end

    return theta_origin, phi_origin, dtheta_dxi1, dphi_dxi2
end

# Initialize auxiliary variables for covariant equations on a torus in Cartesian 3D
function TrixiAtmo.init_auxiliary_node_variables!(auxiliary_variables,
                                                  mesh::P4estMesh{2, 3},
                                                  equations::TrixiAtmo.AbstractCovariantEquations{2,
                                                                                                  3},
                                                  dg,
                                                  elements,
                                                  metric_terms::MetricTermsCovariantTorus,
                                                  bottom_topography)
    @assert equations.global_coordinate_system isa TrixiAtmo.GlobalCartesianCoordinates

    aux_type = eltype(auxiliary_variables.aux_node_vars)
    one_aux = one(aux_type)

    major_radius = convert(aux_type, metric_terms.major_radius)
    minor_radius = convert(aux_type, metric_terms.minor_radius)
    @assert major_radius > zero(aux_type)
    @assert minor_radius > zero(aux_type)

    theta_origin, phi_origin, dtheta_dxi1, dphi_dxi2 = calc_torus_element_map_parameters(mesh,
                                                                                           metric_terms,
                                                                                           aux_type)

    # Exact element-local torus map X(ξ₁, ξ₂) using p4est tree/quadrant metadata
    surface_map_for_element = let theta_origin = theta_origin,
                                  phi_origin = phi_origin,
                                  dtheta_dxi1 = dtheta_dxi1,
                                  dphi_dxi2 = dphi_dxi2,
                                  major_radius = major_radius,
                                  minor_radius = minor_radius,
                                  one_aux = one_aux
        element -> begin
            theta0 = theta_origin[element]
            phi0 = phi_origin[element]
            dtheta = dtheta_dxi1[element]
            dphi = dphi_dxi2[element]

            return (xi1, xi2) -> begin
                theta = theta0 + dtheta * (xi1 + one_aux)
                phi = phi0 + dphi * (xi2 + one_aux)
                torus2cartesian(theta, phi, major_radius, minor_radius)
            end
        end
    end

    init_auxiliary_node_variables_from_map!(auxiliary_variables, elements, mesh,
                                            equations, dg, surface_map_for_element,
                                            bottom_topography)

    return nothing
end

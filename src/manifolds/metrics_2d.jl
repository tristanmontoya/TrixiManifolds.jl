# Shared helpers for metric computations on 2D manifolds embedded in 3D

# Each p4est level halves both reference directions, so the quadrant scale is a scalar
@inline function p4est_quadrant_reference_scale(level::Integer,
                                                ::Type{RealT}) where {RealT <: Real}
    return ldexp(one(RealT), -Int(level))
end

# Generic autodiff helper: compute covariant basis vectors from an arbitrary surface map
# X(ξ₁, ξ₂)
@inline function calc_basis_covariant_autodiff(surface_map, xi1, xi2)
    # Differentiate X(ξ₁, ξ₂) along each reference direction to get tangent vectors
    dxdxi1 = TrixiAtmo.derivative(x -> surface_map(x, xi2), xi1)
    dxdxi2 = TrixiAtmo.derivative(x -> surface_map(xi1, x), xi2)
    return SMatrix{3, 2}(dxdxi1[1], dxdxi1[2], dxdxi1[3],
                         dxdxi2[1], dxdxi2[2], dxdxi2[3])
end

# Compute the covariant metric tensor G_cov = basis_covariantᵀ ⋅ basis_covariant from the
# surface map
@inline function calc_metric_covariant_autodiff(surface_map, xi1, xi2)
    basis_covariant = calc_basis_covariant_autodiff(surface_map, xi1, xi2)
    Gcov = basis_covariant' * basis_covariant
    return SVector(Gcov[1, 1], Gcov[1, 2], Gcov[2, 2])
end

# Compute derivatives of the covariant metric tensor with respect to reference coordinates for Christoffel symbol calculations
@inline function calc_metric_derivatives_autodiff(surface_map, xi1, xi2)
    dGdxi1 = TrixiAtmo.derivative(x -> calc_metric_covariant_autodiff(surface_map, x, xi2),
                                  xi1)
    dGdxi2 = TrixiAtmo.derivative(x -> calc_metric_covariant_autodiff(surface_map, xi1, x),
                                  xi2)
    return dGdxi1, dGdxi2
end

# Use the metric derivatives to compute Christoffel symbols of the second kind for
# covariant equations on a 2D manifold
@inline function calc_christoffel_symbols_covariant(dGdxi1, dGdxi2, Gcon)
    dG11_dxi1, dG12_dxi1, dG22_dxi1 = dGdxi1
    dG11_dxi2, dG12_dxi2, dG22_dxi2 = dGdxi2
    half = convert(eltype(dGdxi1), 0.5)

    # Christoffel symbols of the first kind for 2D metric data
    Gamma_1 = SMatrix{2, 2}(half * dG11_dxi1,
                            half * dG11_dxi2,
                            half * dG11_dxi2,
                            dG12_dxi2 - half * dG22_dxi1)
    Gamma_2 = SMatrix{2, 2}(dG12_dxi1 - half * dG11_dxi2,
                            half * dG22_dxi1,
                            half * dG22_dxi1,
                            half * dG22_dxi2)

    # Raise the first index with G_con to get the symbols of the second kind
    return SVector(Gcon[1, 1] * Gamma_1[1, 1] + Gcon[1, 2] * Gamma_2[1, 1],
                   Gcon[1, 1] * Gamma_1[1, 2] + Gcon[1, 2] * Gamma_2[1, 2],
                   Gcon[1, 1] * Gamma_1[2, 2] + Gcon[1, 2] * Gamma_2[2, 2],
                   Gcon[2, 1] * Gamma_1[1, 1] + Gcon[2, 2] * Gamma_2[1, 1],
                   Gcon[2, 1] * Gamma_1[1, 2] + Gcon[2, 2] * Gamma_2[1, 2],
                   Gcon[2, 1] * Gamma_1[2, 2] + Gcon[2, 2] * Gamma_2[2, 2])
end

# Initialize covariant auxiliary variables from element-local surface maps X(ξ₁, ξ₂)
# `surface_map_for_element(element)` must return a callable `(ξ₁, ξ₂) → SVector{3}`
function init_auxiliary_node_variables_from_map!(auxiliary_variables,
                                                 elements,
                                                 mesh::P4estMesh{2, 3},
                                                 equations::TrixiAtmo.AbstractCovariantEquations{2,
                                                                                                 3},
                                                 dg,
                                                 surface_map_for_element,
                                                 bottom_topography)
    (; node_coordinates) = elements
    (; aux_node_vars) = auxiliary_variables

    zero_aux = zero(eltype(aux_node_vars))
    n_aux = TrixiAtmo.n_aux_node_vars(equations)

    # The minimum number of auxiliary variables for 2D covariant equations is 26, but 
    # more can be allocated for user-defined purposes (e.g. extra physics variables)
    @assert n_aux >= 26

    # Auxiliary variable layout for covariant equations:
    # 1:6 covariant basis, 7:12 contravariant basis, 13 √det(G_cov),
    # 14:16 G_cov, 17:19 G_con, 20 bottom topography, 21:26 Christoffel symbols
    Trixi.@threaded for element in 1:Trixi.ncells(mesh)
        # Define exact element-local surface map from reference element into R³
        surface_map = surface_map_for_element(element)

        for j in eachnode(dg), i in eachnode(dg)
            # Extract reference node positions
            xi1 = dg.basis.nodes[i]
            xi2 = dg.basis.nodes[j]

            # Map to global Cartesian coordinates on the torus for this node/element
            x_node = surface_map(xi1, xi2)

            # Overwrite interpolated element node coordinates with exact map points so all
            # geometry use (metrics, fluxes, callbacks) is evaluated at the same nodes
            node_coordinates[:, i, j, element] .= x_node

            # Covariant/contravariant bases and metric tensors
            basis_covariant = calc_basis_covariant_autodiff(surface_map, xi1, xi2)
            Gcov = basis_covariant' * basis_covariant
            Gcon = inv(Gcov)
            basis_contravariant = Gcon * basis_covariant'

            # Metric derivatives for Christoffel symbols
            dGdxi1, dGdxi2 = calc_metric_derivatives_autodiff(surface_map, xi1, xi2)

            # Fill auxiliary variables for this node/element with exact map data
            aux_node_vars[1:6, i, j, element] = SVector(basis_covariant)
            aux_node_vars[7:12, i, j, element] = SVector(basis_contravariant)
            aux_node_vars[13, i, j, element] = sqrt(det(Gcov))
            aux_node_vars[14:16, i, j, element] = SVector(Gcov[1, 1], Gcov[1, 2],
                                                          Gcov[2, 2])
            aux_node_vars[17:19, i, j, element] = SVector(Gcon[1, 1], Gcon[1, 2],
                                                          Gcon[2, 2])
            aux_node_vars[20, i, j, element] = isnothing(bottom_topography) ?
                                               zero_aux : bottom_topography(x_node)
            aux_node_vars[21:26, i, j, element] = calc_christoffel_symbols_covariant(dGdxi1,
                                                                                     dGdxi2,
                                                                                     Gcon)
            # Fill any remaining auxiliary variables with zeros
            if n_aux > 26
                aux_node_vars[27:n_aux, i, j, element] .= zero_aux
            end
        end
    end

    return nothing
end

@doc raw"""
    P4estMeshTorus2D(trees_per_major_angle, trees_per_minor_angle,
                     major_radius, minor_radius;
                     polydeg, RealT=Float64,
                     initial_refinement_level=0,
                     unsaved_changes=true,
                     p4est_partition_allow_for_coarsening=true)

Create a periodic `P4estMesh{2,3}` representing a torus surface embedded in 3D.
The two logical mesh directions correspond to:
- major angle `theta`: rotation about the global `z` axis
- minor angle `phi`: rotation around the tube cross section

Coordinate convention used by `torus2cartesian`:
- `theta = 0` is on the positive `x` axis and increases toward `+y`.
- `phi = 0` points from the tube centerline outward, away from the `z` axis.
- At `(theta, phi) = (0, 0)`, the point is `(major_radius + minor_radius, 0, 0)`.
- Thus both angles are "aligned with +x" at zero only at that reference location.

Element reference coordinates ``(\xi_1, \xi_2)`` are mapped directly to torus angles
``(\theta, \phi)`` and then to Cartesian coordinates via
```math
\begin{aligned}
\theta &= \frac{2\pi}{N_\theta}\left(k_\theta - 1 + \frac{\xi_1 + 1}{2}\right),\\
\phi &= \frac{2\pi}{N_\phi}\left(k_\phi - 1 + \frac{\xi_2 + 1}{2}\right),\\
x &= (R + r\cos\phi)\cos\theta,\\
y &= (R + r\cos\phi)\sin\theta,\\
z &= r\sin\phi.
\end{aligned}
```
with major radius ``R``, minor radius ``r``, element indices
``k_\theta = 1,\dots,N_\theta``, ``k_\phi = 1,\dots,N_\phi``, and local reference
coordinates ``(\xi_1,\xi_2) \in [-1,1]^2``.
"""
function P4estMeshTorus2D(trees_per_major_angle::Integer,
                          trees_per_minor_angle::Integer,
                          major_radius,
                          minor_radius;
                          polydeg,
                          RealT = Float64,
                          initial_refinement_level = 0,
                          unsaved_changes = true,
                          p4est_partition_allow_for_coarsening = true)
    trees_per_dimension = (trees_per_major_angle, trees_per_minor_angle)
    return P4estMeshTorus2D(trees_per_dimension, major_radius, minor_radius;
                            polydeg = polydeg,
                            RealT = RealT,
                            initial_refinement_level = initial_refinement_level,
                            unsaved_changes = unsaved_changes,
                            p4est_partition_allow_for_coarsening = p4est_partition_allow_for_coarsening)
end

function P4estMeshTorus2D(trees_per_dimension::NTuple{2, <:Integer},
                          major_radius,
                          minor_radius;
                          polydeg,
                          RealT = Float64,
                          initial_refinement_level = 0,
                          unsaved_changes = true,
                          p4est_partition_allow_for_coarsening = true)
    major_radius = RealT(major_radius)
    minor_radius = RealT(minor_radius)

    @assert major_radius > zero(RealT)
    @assert minor_radius > zero(RealT)
    @assert major_radius > minor_radius

    basis = LobattoLegendreBasis(RealT, polydeg)
    nodes = basis.nodes

    n_trees = prod(trees_per_dimension)
    tree_node_coordinates = Array{RealT, 4}(undef, 3,
                                            ntuple(_ -> length(nodes), 2)...,
                                            n_trees)

    # Map (element id, local (ξ₁, ξ₂)) → (θ, φ) → (x, y, z)
    calc_torus_tree_node_coordinates!(tree_node_coordinates, nodes, trees_per_dimension,
                                      major_radius, minor_radius)

    periodicity = (true, true)
    connectivity = Trixi.connectivity_structured(trees_per_dimension..., periodicity)
    p4est = Trixi.new_p4est(connectivity, initial_refinement_level)

    boundary_names = fill(Symbol("---"), 2 * 2, n_trees)

    return P4estMesh{2}(p4est, tree_node_coordinates, nodes,
                        boundary_names, "", unsaved_changes,
                        p4est_partition_allow_for_coarsening)
end

# Calculate torus coordinates for each node based on element index and local reference
# coordinates
function calc_torus_tree_node_coordinates!(tree_node_coordinates::Array{RealT, 4},
                                           nodes,
                                           trees_per_dimension::NTuple{2, <:Integer},
                                           major_radius,
                                           minor_radius) where {RealT}
    trees_per_major_angle, trees_per_minor_angle = trees_per_dimension
    linear_indices = LinearIndices(trees_per_dimension)

    # Type conversions for constants used in the mapping formulas
    two_pi = convert(RealT, 2pi)
    half = convert(RealT, 0.5)
    one_rt = one(RealT)

    for k_phi in 1:trees_per_minor_angle, k_theta in 1:trees_per_major_angle
        tree_id = linear_indices[k_theta, k_phi]

        for j in eachindex(nodes), i in eachindex(nodes)
            xi_1 = nodes[i]
            xi_2 = nodes[j]

            # Map to global torus angles on this element, which are affine in local
            # reference coordinates
            theta = two_pi *
                    ((k_theta - 1) + (xi_1 + one_rt) * half) /
                    trees_per_major_angle
            phi = two_pi *
                  ((k_phi - 1) + (xi_2 + one_rt) * half) /
                  trees_per_minor_angle

            # Map from torus angles to Cartesian coordinates for this node/element
            tree_node_coordinates[:, i, j, tree_id] .= torus2cartesian(theta, phi,
                                                                       major_radius,
                                                                       minor_radius)
        end
    end

    return nothing
end

@doc raw"""
    torus2cartesian(theta, phi, major_radius, minor_radius)

Map global torus coordinates ``(\theta, \phi)`` to global Cartesian coordinates
``(x, y, z)`` using
```math
\begin{aligned}
x &= (R + r\cos\phi)\cos\theta,\\
y &= (R + r\cos\phi)\sin\theta,\\
z &= r\sin\phi,
\end{aligned}
```
where `R = major_radius` and `r = minor_radius`.
"""
@inline function torus2cartesian(theta, phi, major_radius, minor_radius)
    s_theta, c_theta = sincos(theta)
    s_phi, c_phi = sincos(phi)

    # Parametrization with symmetry axis along +z:
    # x = (R + r cos(φ)) cos(θ)
    # y = (R + r cos(φ)) sin(θ)
    # z = r sin(φ)
    ring_radius = major_radius + minor_radius * c_phi
    x = ring_radius * c_theta
    y = ring_radius * s_theta
    z = minor_radius * s_phi

    return SVector(x, y, z)
end

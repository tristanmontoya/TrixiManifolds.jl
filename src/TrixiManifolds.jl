module TrixiManifolds

using Trixi
using TrixiAtmo
using LinearAlgebra: det, inv

export examples_dir
export P4estMeshTorus2D, MetricTermsCovariantTorus, torus2cartesian

const examples_dir = joinpath(@__DIR__, "..", "examples")

include("manifolds/manifolds.jl")

end

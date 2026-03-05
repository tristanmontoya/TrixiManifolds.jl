# TrixiManifolds.jl
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

**Note: This repository is still in its alpha stage and anything might change at any time and without warning.**

`TrixiManifolds.jl` provides extensions for solving PDEs on curved manifolds using the [Trixi.jl](https://github.com/trixi-framework/Trixi.jl) ecosystem.

## Included features
- Torus manifold support: custom periodic torus mesh construction and covariant metric-term initialization, demonstrated in `examples/elixir_advection_torus.jl`.
- Sphere support included from [TrixiAtmo.jl](https://github.com/trixi-framework/TrixiAtmo.jl): covariant advection on a cubed sphere demonstrated in `examples/elixir_advection_cubed_sphere.jl`.

## Installation
First, make sure you have [Julia](https://julialang.org/downloads/) installed (the code was tested with Julia v1.12). Then, assuming you're on Linux or MacOS, run the following commands:

```bash
git clone git@github.com:tristanmontoya/TrixiManifolds.jl.git
cd TrixiManifolds.jl
mkdir run
cd run
julia --project=. -e 'using Pkg; Pkg.develop(PackageSpec(path=".."))'
julia --project=. -e 'using Pkg; Pkg.add(["Trixi", "TrixiAtmo", "OrdinaryDiffEq", "Trixi2Vtk"])'
```

## Basic usage
From the `run` directory, start Julia with `julia --project=.` and run one of the manifold advection examples. Sphere support is included out of the box:

```julia
julia> using Trixi
julia> trixi_include("../examples/elixir_advection_cubed_sphere.jl")
```

or

```julia
julia> trixi_include("../examples/elixir_advection_torus.jl")
```

These examples are configured for `T = 4` with output interval `0.02`. They already convert output to `.vtu` at the end. If needed, you can run:

```julia
julia> using Trixi2Vtk
julia> trixi2vtk("out/advection_torus/solution_*.h5", output_directory="out/advection_torus")
julia> trixi2vtk("out/advection_cubed_sphere/solution_*.h5", output_directory="out/advection_cubed_sphere")
```

## Headless animations (no GUI)
Use the script in `scripts/` to render animations directly from terminal. It auto-reexecutes via `pvbatch` when run with normal Python.

Prerequisites:
- Python 3
- ParaView installed with `pvbatch` available (either on `PATH` or via `PVBATCH=/path/to/pvbatch`)
- `ffmpeg` on `PATH` for movie output (`.mp4`, `.avi`, etc.); PNG plot output does not require `ffmpeg`
- VTU files in `out/` (for example from `trixi2vtk("out/solution_*.h5", output_directory="out/")`)

```bash
python scripts/render_paraview_animation.py \
  --input "out/advection_torus/solution_*.vtu" \
  --output out/advection_torus/animation.mp4 \
  --field h
```

```bash
python scripts/render_paraview_animation.py \
  --input "out/advection_cubed_sphere/solution_*.vtu" \
  --output out/advection_cubed_sphere/animation.mp4 \
  --field h \
  --camera isometric
```

## Headless ParaView plots (PNG)
To render a single still plot, set `--start-frame` and `--end-frame` to the same index and use a quoted wildcard output like `plot_*.png`. The script normalizes this to ParaView's filename seed and writes files matching `plot_*.png` (for example, `plot_.0100.png`).

```bash
python scripts/render_paraview_animation.py \
  --input "out/advection_torus/solution_*.vtu" \
  --output "out/advection_torus/plot_*.png" \
  --field h \
  --start-frame 100 \
  --end-frame 100
```

```bash
python scripts/render_paraview_animation.py \
  --input "out/advection_cubed_sphere/solution_*.vtu" \
  --output "out/advection_cubed_sphere/plot_*.png" \
  --field h \
  --start-frame 100 \
  --end-frame 100
```

## License
This code is released under the [MIT license](https://opensource.org/license/mit).

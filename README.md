
# 💧 FluidProject — 2D Stable Fluids

A real-time interactive 2D fluid simulation built on **Jos Stam's Stable Fluids** algorithm (SIGGRAPH 1999). The solver is unconditionally stable, allowing large timesteps while producing physically plausible smoke/dye advection, velocity diffusion, and pressure projection — all running in real time.

***

## 📸 Demo

> _Add a GIF or screenshot of the simulation here._
> Example: `![Fluid Simulation Demo](demo.gif)`

***

## 🧠 Algorithm Overview

The simulation solves the **incompressible Navier-Stokes equations** on a uniform 2D grid using Stam's operator-splitting approach:

Each timestep runs four stages:

1. **Add Sources** — inject velocity or dye density from user input
2. **Diffuse** — spread velocity/density using an implicit (Gauss-Seidel) solver
3. **Advect** — semi-Lagrangian backtracing: for each cell, trace backward along the velocity field and interpolate
4. **Project** — enforce incompressibility (divergence-free velocity field) via pressure solve

The key insight from Stam's paper is the **semi-Lagrangian advection** scheme, which is unconditionally stable for any timestep size — eliminating the CFL constraint that plagues explicit methods.

### Core Math

The velocity field **u** evolves by:

$$\frac{\partial \mathbf{u}}{\partial t} = -(\mathbf{u} \cdot \nabla)\mathbf{u} + \nu \nabla^2 \mathbf{u} + \mathbf{f}$$

Where:
- `ν` = kinematic viscosity
- `f` = external forces (mouse input)
- The divergence-free constraint: `∇ · u = 0`

***

## 🚀 Getting Started

### Prerequisites

- C++ compiler (g++ / clang++)
- OpenGL + GLFW + GLEW (or equivalent window/rendering lib)
- CMake (optional, if using CMakeLists.txt)

### Build & Run

```bash
git clone https://github.com/VeliuSami/FluidProject-2DStableFluids.git
cd FluidProject-2DStableFluids
mkdir build && cd build
cmake ..
make
./FluidProject
```

Or compile manually:
```bash
g++ -std=c++17 -o FluidProject main.cpp -lGL -lGLFW -lGLEW
./FluidProject
```

***

## 🎮 Controls

| Input | Action |
|-------|--------|
| `Left Click + Drag` | Add velocity force |
| `Right Click + Drag` | Add dye / density |
| `R` | Reset simulation |
| `V` | Toggle velocity field visualization |
| `D` | Toggle density visualization |
| `+` / `-` | Increase / decrease viscosity |
| `ESC` | Quit |

> _Update these keybindings to match your actual implementation._

***

## ⚙️ Parameters

These can be tweaked in `config.h` or at the top of `main.cpp`:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `N` | Grid resolution (N×N) | `128` |
| `dt` | Timestep | `0.1` |
| `diff` | Diffusion rate | `0.0` |
| `visc` | Viscosity | `0.0` |
| `iter` | Gauss-Seidel iterations | `20` |

***

## 🙈 Hiding Files from Git

To hide (exclude) any file from being tracked by Git, add it to `.gitignore`.

### Step 1 — Create or edit `.gitignore`

```bash
# In the root of your repo
touch .gitignore
```

### Step 2 — Add patterns for files you want hidden

Open `.gitignore` in any editor and add:

```gitignore
# Ignore a specific file
secrets.txt
config/api_keys.h

# Ignore all files with a specific extension
*.log
*.o
*.out
*.exe

# Ignore an entire folder
build/
__pycache__/
.vscode/

# Ignore everything EXCEPT a specific file (use ! to negate)
*.env
!.env.example
```

### Step 3 — Apply to already-tracked files

If a file is **already committed** and you want Git to stop tracking it:

```bash
git rm --cached filename.txt        # stop tracking one file
git rm --cached -r build/           # stop tracking a whole folder
git commit -m "Remove tracked files now in .gitignore"
```

> ⚠️ `git rm --cached` does NOT delete the file locally — it only tells Git to stop watching it.

### Step 4 — Verify it's hidden

```bash
git status          # file should no longer appear as modified/untracked
git check-ignore -v filename.txt    # debug: confirm the rule that's hiding it
```

### Common Patterns for This Project

```gitignore
# Build artifacts
build/
*.o
*.out
*.exe
FluidProject

# IDE files
.vscode/
.idea/
*.xcodeproj

# OS files
.DS_Store
Thumbs.db

# Sensitive config
config/secrets.h
*.env
```

***

## 📚 References

- Jos Stam, [*Stable Fluids*](https://pages.cs.wisc.edu/~chaol/data/cs777/stam-stable_fluids.pdf), SIGGRAPH 1999
- Mark Harris, [*Fast Fluid Dynamics Simulation on the GPU*](https://developer.nvidia.com/gpugems/gpugems/part-vi-beyond-triangles/chapter-38-fast-fluid-dynamics-simulation-gpu), GPU Gems
- Robert Bridson, [*Fluid Simulation for Computer Graphics*](https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf)

***

## 📄 License

MIT License — feel free to use, modify, and distribute.

***

_Built by [@VeliuSami](https://github.com/VeliuSami)_

### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 180d13c8-ffcf-11ee-3ce4-75582cab0f99
import Pkg; Pkg.activate("..")

# ╔═╡ 6a729112-8817-49ff-955d-f8393f9bc375
using PlutoUI, Plots, LinearAlgebra

# ╔═╡ 0a49b2e8-a425-46b2-a392-1a3a55e7198f
# use functions from last tutorials (4,5)

# mom2Vel
# calcFluxConCDS
# calcFluxDif
# difMom2D
# poissonSolver2D
# corrMom2D

# ╔═╡ daf885c9-6e69-4e78-9bdf-def4e6f4e766
# Passive scalar transport in 2D

# 1. prescribe a momentum field and a constant density
# 2. compute velocity at faces (mom2Vel)
# 3. compute the convective and diffusive fluxes (calcFluxConCds, calcFluxDif) at faces
# 4. sum the computed fluxes up at the cell center
# 5. advance the solution in time by adding the flux at the cell center to the previous solution
# 6. repeat last 3 steps until the desired number of timesteps is reached (as 1D code in lec3)
# 7. testcase: jet

# ╔═╡ 08667b2d-563e-4615-9598-3e9da7c6121e
# Simple flow solver in 2D

# 1. extend "passive scalar transport" to also consider the transport of momentum in x- and y-direction
# 2. after advanced the solution in time, you need to do the pressure correction in addition for the momentum (divMom2D, poissonSolver2D, corrMom2D)
# 3. at the end of the time step, you need to update the velocity at faces (mom2Vel)
# 4. try effect of different parameterse, discretisation schemes for convection, reynolds numbers, grid resolution

# ╔═╡ Cell order:
# ╠═180d13c8-ffcf-11ee-3ce4-75582cab0f99
# ╠═6a729112-8817-49ff-955d-f8393f9bc375
# ╠═0a49b2e8-a425-46b2-a392-1a3a55e7198f
# ╠═daf885c9-6e69-4e78-9bdf-def4e6f4e766
# ╠═08667b2d-563e-4615-9598-3e9da7c6121e

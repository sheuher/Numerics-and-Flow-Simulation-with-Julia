### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ ccc43210-ffcd-11ee-0ae6-9138d9dd2fb8
import Pkg; Pkg.activate("..")

# ╔═╡ dcb898a7-fe46-47bf-ae4c-4e7584573edc
using PlutoUI, Plots, LinearAlgebra

# ╔═╡ c9901a62-25b8-4fba-a806-d9d40ad051cd
# 1. Projection method in 1D

# ╔═╡ 6b8f925b-cd26-487a-ae16-c17e49db890b
function divMom1D()
	# calc resulting divergence field for a given momentum prediction field using CDS
end

# ╔═╡ 1b209e76-af07-449e-bb40-2b8cd7225728
function poissonSolver1D()
	# solve poisson eq. with iterative jacobi method for given divergence field, apply zero gradient BC at the inlet and fix pressure = 0 at outlet
end

# ╔═╡ 18abc503-8936-412f-ad24-395b9c4b44d4
function corrMom1D()
	# corrects momentum prediction field based on given pressure field
	# use CDS for pressure gradient
end

# ╔═╡ a1e1d09b-0528-4586-a035-9eae76b6a557
# program that perform the 3 former steps to correct a given momentum prediction field

# ╔═╡ 05ce9937-5f4e-405d-b27e-0d8c3707c743
# test against lecture

# ╔═╡ 5f2d47bb-6104-4572-9290-a88111173886
# use analytic solutions to validate program

# ╔═╡ b3c39e35-ce7d-42b2-8e54-5e53d34bed2f
### use ghostcells for bc and stick to the naming conventions

# ╔═╡ 2bda1f1e-eb89-4698-86b4-677d4c1bdcd9


# ╔═╡ dfbe2575-9de8-4711-9073-28745d6b5d19


# ╔═╡ d4d20e5c-9809-466e-a4fd-8aab5001a309
# 2. Projection method in 2D

# ╔═╡ 744664e2-eee8-466e-b1ec-8a1d0b7533bc
## extend carefully from 1D to 2D

# ╔═╡ 265c8c4f-eac4-49c5-ad03-86fbe03321a0
## BCs for pressure = 0 gradient at inlet; fixed pressure = 0 at all 3

# ╔═╡ efdf311d-29c2-4e2e-99a3-a623c21268eb
## test: compare pressure field from 2D example 5 to the sum of 2D examples 1-4

# ╔═╡ Cell order:
# ╠═ccc43210-ffcd-11ee-0ae6-9138d9dd2fb8
# ╠═dcb898a7-fe46-47bf-ae4c-4e7584573edc
# ╠═c9901a62-25b8-4fba-a806-d9d40ad051cd
# ╠═6b8f925b-cd26-487a-ae16-c17e49db890b
# ╠═1b209e76-af07-449e-bb40-2b8cd7225728
# ╠═18abc503-8936-412f-ad24-395b9c4b44d4
# ╠═a1e1d09b-0528-4586-a035-9eae76b6a557
# ╠═05ce9937-5f4e-405d-b27e-0d8c3707c743
# ╠═5f2d47bb-6104-4572-9290-a88111173886
# ╠═b3c39e35-ce7d-42b2-8e54-5e53d34bed2f
# ╠═2bda1f1e-eb89-4698-86b4-677d4c1bdcd9
# ╠═dfbe2575-9de8-4711-9073-28745d6b5d19
# ╠═d4d20e5c-9809-466e-a4fd-8aab5001a309
# ╠═744664e2-eee8-466e-b1ec-8a1d0b7533bc
# ╠═265c8c4f-eac4-49c5-ad03-86fbe03321a0
# ╠═efdf311d-29c2-4e2e-99a3-a623c21268eb

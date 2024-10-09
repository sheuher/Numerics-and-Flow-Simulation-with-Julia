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
function divMom1D!(divPred, rhouP, delta)
	# calc resulting divergence field for a given momentum prediction field using CDS

	nG = 1
	ImaAll = length(divPred)
	Ifi = nG +1; Ifim = Ifi-1; Ifip = Ifi +1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1

	divPred[Ifi:Ila] .= 0.5 /delta *(rhouP[Ifip:Ilap] .- rhouP[Ifim:Ilam]) 

	return nothing
end

# ╔═╡ 1b209e76-af07-449e-bb40-2b8cd7225728
function poissonSolver1D!(p, divPred, dt, delta, nItMax, epsMax)
	# solve poisson eq. with iterative jacobi method for given divergence field, apply zero gradient BC at the inlet and fix pressure = 0 at outlet

	nG = 1
	ImaAll = length(divPred)
	Ifi = nG +1; Ifim = Ifi-1; Ifip = Ifi +1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1

	RHS = divPred /dt
	nIt = 0; eps = 10.0
	while ( (nIt <nItMax) && (eps>epsMax) )
		pold = copy(p)
		pold[Ifim] = 101
		pold[Ilap] = 0
		p[Ifi:Ila] = 0.5 *( pold[Ifim:Ilam] +pold[Ifip:Ilap] -RHS[Ifi:Ila] *delta^2 )
		p[Ifim] = p[Ifi]
		p[Ilap] = -p[Ila]
		eps = maximum( abs.(p[Ifi:Ila] -pold[Ifi:Ila]) )
		nIt = nIt + 1
	end
	
	return (nIt=nIt, eps=eps)
end

# ╔═╡ 18abc503-8936-412f-ad24-395b9c4b44d4
function corrMom1D!(rhou, rhouP, p, dt, delta)
	# corrects momentum prediction field based on given pressure field
	# use CDS for pressure gradient

	nG = 1
	ImaAll = length(rhou)
	Ifi = nG +1; Ifim = Ifi-1; Ifip = Ifi +1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1

	rhou[Ifi:Ila] = rhouP[Ifi:Ila] .- dt * 0.5 /delta *(p[Ifip:Ilap] .-p[Ifim:Ilam]) 

	return nothing
end

# ╔═╡ 74f5b067-a79c-4d69-806b-ef9758dc68b9
let
	x = 3π; CFL = 0.25; U = 0.25; D = 1e-5

	nItMax = 100000; epsMax = 1e-9
	
	Ima = 100; nG = 1
	Ifim = 1; Ifi = nG +1; 
	Ila = Ima +nG; Ilap = Ila + 1
	
	rhouP = Vector{Float64}(undef, Ilap)
	divPred = zeros(Ilap); p = zeros(Ilap)
	rhou = zeros(Ilap)

	println("x = $x, U = $U, CFL = $(CFL), D = $(D)")
	println("Ima = $(Ima)")

	DeltaX = x / Ima
	DeltaT = CFL *DeltaX /U
	println("DeltaT = $(DeltaT)")
	println("DeltaX = $(DeltaX)")
	println("tStep = $(step)")
	
	rhouP[Ifi:Ila] .= LinRange(DeltaX/2,x - DeltaX/2,Ima) .|> t->cos(t)
	rhouP[Ifim] = rhouP[Ila]; rhouP[Ilap] = rhouP[Ifi]

	# pressure initialisation
	p[Ifim] = 101
	p[Ilap] = 0


	
	divMom1D!(divPred, rhouP, DeltaX)
	poissonError = poissonSolver1D!(p, divPred, DeltaT, DeltaX, nItMax, epsMax); println(poissonError)
	corrMom1D!(rhou, rhouP, p, DeltaT, DeltaX)


	
	p1 = plot(rhou[Ifi:Ila], label="corrected momentum")
	plot!(p1, - DeltaT * 0.5 /DeltaX *(p[Ifi+1:Ilap] .-p[Ifim:Ila-1]), label="correction")
	plot!(p1, rhouP[Ifi:Ila], label="predicted momentum")
	p2 = plot(p, label="jacobi pressure at iter steps = $(poissonError[1])")
	plot(p1,p2,layout=(2,1),size=(750,900))
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

# ╔═╡ 172fbd89-20b9-4eed-b5fc-d40d2a86f01e


# ╔═╡ 3456f4e4-32f0-42ea-a9a7-b57192f1c24f


# ╔═╡ eaef9a1d-bcd3-45d4-9d29-bc7cbecba568
"""
	divMom2D!(divPred, rhouP, delta)
	#heatmap(divPred)
	poissonSolver2D!(p, divPred, dt, delta, 10000, epsMax) |> println
	heatmap(divPred[Ifi:Ila, Jfi:Jla])
	heatmap(p[Ifi:Ila, Jfi:Jla])
	"""

# ╔═╡ aff8d611-92c7-43de-b559-ed672901a33e
function divMom2D!(divPred, rhouP, delta)
	# calc resulting divergence field for a given momentum prediction field using CDS

	nG = 1
	ImaAll, JmaAll = size(divPred)
	Ifi = nG +1; Ifim = Ifi-1; Ifip = Ifi +1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG +1; Jfim = Jfi-1; Jfip = Jfi +1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1

	rhouPE = @view rhouP[Ifip:Ilap, Jfi:Jla, 1]
	rhouPW = @view rhouP[Ifim:Ilam, Jfi:Jla, 1]
	rhouPN = @view rhouP[Ifi:Ila, Jfip:Jlap, 2]
	rhouPS = @view rhouP[Ifi:Ila, Jfim:Jlam, 2]

	divPred[Ifi:Ila, Jfi:Jla] .= 0.5 /delta *( rhouPE .- rhouPW .+ rhouPN .- rhouPS ) 

	return nothing
end

# ╔═╡ 7f6ca118-ae7c-4790-b5b0-b2e7a05a6b49
function poissonSolver2D!(p, divPred, dt, delta, nItMax, epsMax)
	# solve poisson eq. with iterative jacobi method for given divergence field, apply zero gradient BC at the inlet and fix pressure = 0 at outlet

	nG = 1
	ImaAll, JmaAll = size(divPred)
	Ifi = nG +1; Ifim = Ifi-1; Ifip = Ifi +1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG +1; Jfim = Jfi-1; Jfip = Jfi +1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1
	
	RHS = divPred /dt
	nIt = 0; eps = 10.0
	### while ( (nIt <nItMax) && (eps>epsMax) )
	for _ in 1:nItMax
		pold = copy(p)
		pE = @view pold[Ifip:Ilap, Jfi:Jla]
		pW = @view pold[Ifim:Ilam, Jfi:Jla]
		pN = @view pold[Ifi:Ila, Jfip:Jlap]
		pS = @view pold[Ifi:Ila, Jfim:Jlam]
		p[Ifi:Ila, Jfi:Jla] = 0.25 *(  pE + pW + pN + pS -RHS[Ifi:Ila, Jfi:Jla]*delta^2 )
		p[Ifim, Jfi:Jla] .= p[Ifi, Jfi:Jla]
		eps = maximum( abs.(p[Ifi:Ila] -pold[Ifi:Ila]) )
		#if eps < epsMax
		#	break
		#end
		nIt = nIt + 1
	end
	return (nIt=nIt, eps=eps)
end

# ╔═╡ 1128ed32-e189-4dd4-993c-8b71b717c695
# EXAMPLE 1

let
	Ima = 50; Jma = 50; delta = 0.1
	dt = .1; epsMax = 1e-6

	nG = 1
	ImaAll = Ima+2nG; JmaAll = Jma+2nG

	rhouP = zeros(ImaAll, JmaAll, 2)
	divPred = zeros(ImaAll, JmaAll)
	p = ones(ImaAll, JmaAll)
	
	rhouP[:,:,1] .= 1.0
	rhouP[:,:,2] .= 0
	
	divMom2D!(divPred, rhouP, delta)
	poissonSolver2D!(p, divPred, dt, delta, 10000, epsMax) |> println
	heatmap(p)
	
end

# ╔═╡ 135259e0-7754-4b1c-97e5-e8d3140bcdd8
# EXAMPLE 2

let
	Ima = 50; Jma = 50; delta = 0.1
	dt = .1; epsMax = 1e-6

	nG = 1
	ImaAll = Ima+2nG; JmaAll = Jma+2nG
	Ifi = nG +1; Ifim = Ifi-1; Ifip = Ifi +1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG +1; Jfim = Jfi-1; Jfip = Jfi +1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1

	rhouP = zeros(ImaAll, JmaAll, 2)
	x = zeros(ImaAll, JmaAll, 2)
	divPred = zeros(ImaAll, JmaAll)
	p = ones(ImaAll, JmaAll)

	tmpx = LinRange(-delta, delta*Ima + delta, ImaAll)
	tmpy = LinRange(-delta, delta*Jma + delta, JmaAll)
	x[:,:,1] = tmpx * ones(length(tmpy))'
	x[:,:,2] = ones(length(tmpx)) * tmpy'

	rhouP[:,:,1] = x[:,:,2] .- (Jma/2 +1) *delta
	rhouP[:,:,2] = -x[:,:,1] .+ (Ima/2 +1) *delta 

	#quiver(x[:,:,1], x[:,:,2], quiver=(rhouP[:,:,1]*1e-3, rhouP[:,:,2]*1e-3))

	#BC zero gradient at Ifi
	rhouP[Ifim,:,1] .= rhouP[Ifi,:,1]
	rhouP[Ifim,:,2] .= rhouP[Ifi,:,2]
	rhouP
	
	divMom2D!(divPred, rhouP, delta)
	#heatmap(divPred)
	poissonSolver2D!(p, divPred, dt, delta, 10000, epsMax) |> println
	heatmap(divPred[Ifi:Ila, Jfi:Jla])
	heatmap(p[Ifi:Ila, Jfi:Jla])
	
end

# ╔═╡ ab0d9908-d262-4fab-b12b-43c0e063ff3f
# EXAMPLE 3

let
	Ima = 50; Jma = 50; delta = 0.1
	dt = .1; epsMax = 1e-6

	nG = 1
	ImaAll = Ima+2nG; JmaAll = Jma+2nG
	Ifi = nG +1; Ifim = Ifi-1; Ifip = Ifi +1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG +1; Jfim = Jfi-1; Jfip = Jfi +1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1

	rhouP = zeros(ImaAll, JmaAll, 2)
	x = zeros(ImaAll, JmaAll, 2)
	divPred = zeros(ImaAll, JmaAll)
	p = ones(ImaAll, JmaAll)

	tmpx = LinRange(-delta, delta*Ima + delta, ImaAll)
	tmpy = LinRange(-delta, delta*Jma + delta, JmaAll)
	x[:,:,1] = tmpx * ones(length(tmpy))'
	x[:,:,2] = ones(length(tmpx)) * tmpy'

	#rhouP initialisation
	rhouP[:,:,1] .= 1. .- x[:,:,1] /(Ima*delta)
	rhouP[:,:,2] .= 0.
	 

	#p initialisation
	#p[Ifim,:] .= 101
	
	quiver(x[:,:,1], x[:,:,2], quiver=(rhouP[:,:,1], rhouP[:,:,2]))

	#BC zero gradient at Ifi
	
	divMom2D!(divPred, rhouP, delta)
	#divPred[Ifim,:,] .= divPred[Ifi,:,]
	poissonSolver2D!(p, divPred, dt, delta, 10000, epsMax) |> println
	heatmap(divPred[Ifi:Ila, Jfi:Jla], size=(750,500))
	heatmap(p[Ifi:Ila, Jfi:Jla])
end

# ╔═╡ 736a5b64-ca9d-4ca6-80d7-b6154df2b68e
# EXAMPLE 4

let
	Ima = 50; Jma = 50; delta = 0.1
	dt = .1; epsMax = 1e-6

	nG = 1
	ImaAll = Ima+2nG; JmaAll = Jma+2nG
	Ifi = nG +1; Ifim = Ifi-1; Ifip = Ifi +1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG +1; Jfim = Jfi-1; Jfip = Jfi +1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1

	rhouP = zeros(ImaAll, JmaAll, 2)
	x = zeros(ImaAll, JmaAll, 2)
	divPred = zeros(ImaAll, JmaAll)
	p = ones(ImaAll, JmaAll)

	tmpx = LinRange(-delta, delta*Ima + delta, ImaAll)
	tmpy = LinRange(-delta, delta*Jma + delta, JmaAll)
	x[:,:,1] = tmpx * ones(length(tmpy))'
	x[:,:,2] = ones(length(tmpx)) * tmpy'

	#rhouP initialisation
	rhouP[:,:,:] .= 0.
	is=convert(Int64,Ima/2-5):convert(Int64,Ima/2+5)
	js=convert(Int64,Jma/2-5):convert(Int64,Jma/2+5)
	rhouP[is,js, 2] .= 1.
	
	#p initialisation
	#p[Ifim,:] .= 101
	
	quiver(x[:,:,1], x[:,:,2], quiver=(rhouP[:,:,2], rhouP[:,:,1]))

	#BC zero gradient at Ifi
	
	divMom2D!(divPred, rhouP, delta)
	poissonSolver2D!(p, divPred, dt, delta, 10000, epsMax) |> println
	heatmap(divPred[Ifi:Ila, Jfi:Jla], size=(750,500))
	heatmap(p[Ifi:Ila, Jfi:Jla])
end

# ╔═╡ 31ecf6a0-c5cf-439e-accb-e702d73d16e5
function corrMom2D!(rhou, rhouP, p, dt, delta)
	# corrects momentum prediction field based on given pressure field
	# use CDS for pressure gradient

	nG = 1
	ImaAll, JmaAll = size(rhou)
	Ifi = nG +1; Ifim = Ifi-1; Ifip = Ifi +1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG +1; Jfim = Jfi-1; Jfip = Jfi +1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1

	rhou[Ifi:Ila, Jfi:Jla, 1] = rhouP[Ifi:Ila, Jfi:Jla, 1] .- dt * 0.5 /delta *(p[Ifip:Ilap, Jfi:Jla] .-p[Ifim:Ilam, Jfi:Jla]) 

	rhou[Ifi:Ila, Jfi:Jla, 2] = rhouP[Ifi:Ila, Jfi:Jla, 2] .- dt * 0.5 /delta *(p[Ifi:Ila, Jfip:Jlap] .-p[Ifi:Ila, Jfim:Jlam]) 
	
	return nothing
end

# ╔═╡ 350980f1-8bf7-4191-ae79-21699c175f46
# EXAMPLE 5

let
	Ima = 50; Jma = 50; delta = 0.1
	dt = .1; epsMax = 1e-6

	nG = 1
	ImaAll = Ima+2nG; JmaAll = Jma+2nG
	Ifi = nG +1; Ifim = Ifi-1; Ifip = Ifi +1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG +1; Jfim = Jfi-1; Jfip = Jfi +1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1

	rhouP = zeros(ImaAll, JmaAll, 2)
	x = zeros(ImaAll, JmaAll, 2)
	divPred = zeros(ImaAll, JmaAll)
	p = ones(ImaAll, JmaAll)
	rhou = zeros(ImaAll, JmaAll, 2)

	tmpx = LinRange(-delta, delta*Ima + delta, ImaAll)
	tmpy = LinRange(-delta, delta*Jma + delta, JmaAll)
	
	x[:,:,1] = tmpx * ones(length(tmpy))'
	x[:,:,2] = ones(length(tmpx)) * tmpy'

	#rhouP initialisation
	rhouP[:,:,1] .= 1.0
	rhouP[:,:,1] .+= x[:,:,2] .- (Jma/2 +1) *delta
	rhouP[:,:,1] .+= 1. .- x[:,:,1] /(Ima*delta)
	
	rhouP[:,:,2] = -x[:,:,1] .+ (Ima/2 +1) *delta
	is=convert(Int64,Ima/2-5):convert(Int64,Ima/2+5)
	js=convert(Int64,Jma/2-5):convert(Int64,Jma/2+5)
	rhouP[is,js, 2] .+= 1.
	
	#quiver(x[:,:,1], x[:,:,2], quiver=(rhouP[:,:,2], rhouP[:,:,1]))

	#BC zero gradient at Ifi
	
	divMom2D!(divPred, rhouP, delta)
	poissonSolver2D!(p, divPred, dt, delta, 10000, epsMax) |> println
	corrMom2D!(rhou, rhouP, p, dt, delta)
	#heatmap(divPred[Ifi:Ila, Jfi:Jla], size=(750,500))
	heatmap(p[Ifi:Ila, Jfi:Jla], clims=(0,8), cmap=:rainbow, size=(750,600))
end

# ╔═╡ cd984d25-d10a-4f3f-be41-28ea10a1c9d8


# ╔═╡ 65ef231b-c3d7-49f0-8b8c-c2992242e2e0
let
	N = 12
	divPred = zeros(N,N)
	p = copy(divPred)

	p[1,:].=100
	p[12,:].=250
	p[:,1].=50
	p[:,12].=300
	
	dt = 0.01
	delta = 0.1
	poissonSolver2D!(p, divPred, 1, 1, 10000, 1e-10) |> println
	
	nG = 1
	ImaAll, JmaAll = size(divPred)
	Ifi = nG +1; Ifim = Ifi-1; Ifip = Ifi +1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG +1; Jfim = Jfi-1; Jfip = Jfi +1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1
	nIt = 0; eps = 10.0
	
	#pNew = copy(p)
	#for _ in 1:1000
	#	pNew[Ifi:Ila, Jfi:Jla] = 0.5 *( pNew[Ifip:Ilap, Jfi:Jla] + pNew[Ifim:Ilam, Jfi:Jla] + pNew[Ifi:Ila, Jfim:Jlam] +  pNew[Ifi:Ila, Jfip:Jlap] )
	#end
	p
end

# ╔═╡ Cell order:
# ╠═ccc43210-ffcd-11ee-0ae6-9138d9dd2fb8
# ╠═dcb898a7-fe46-47bf-ae4c-4e7584573edc
# ╠═c9901a62-25b8-4fba-a806-d9d40ad051cd
# ╠═74f5b067-a79c-4d69-806b-ef9758dc68b9
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
# ╠═172fbd89-20b9-4eed-b5fc-d40d2a86f01e
# ╠═1128ed32-e189-4dd4-993c-8b71b717c695
# ╠═135259e0-7754-4b1c-97e5-e8d3140bcdd8
# ╠═ab0d9908-d262-4fab-b12b-43c0e063ff3f
# ╠═736a5b64-ca9d-4ca6-80d7-b6154df2b68e
# ╠═350980f1-8bf7-4191-ae79-21699c175f46
# ╠═3456f4e4-32f0-42ea-a9a7-b57192f1c24f
# ╠═eaef9a1d-bcd3-45d4-9d29-bc7cbecba568
# ╠═aff8d611-92c7-43de-b559-ed672901a33e
# ╠═7f6ca118-ae7c-4790-b5b0-b2e7a05a6b49
# ╠═31ecf6a0-c5cf-439e-accb-e702d73d16e5
# ╠═cd984d25-d10a-4f3f-be41-28ea10a1c9d8
# ╠═65ef231b-c3d7-49f0-8b8c-c2992242e2e0

### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 7a635be6-ff7c-11ee-2cdd-afdcd8ccc93d
import Pkg; Pkg.activate("..")

# ╔═╡ 718f1af5-c5ea-4b78-b748-eae0c996cf97
using PlutoUI, Plots, LinearAlgebra

# ╔═╡ ad19c163-839b-47af-a528-4df2499a8bd6
function mom2vel!(u, v, rhou, rhov)

	"""
	calculate velocity fields u[:,:] and v[:,:] on the CELL SURFACES
	from the momentum fields rhou[:,:] and rhov[:,:] stored on cell centre
	"""
	nG = 1

	ImaAll, JmaAll = size(u)
	Ifi = nG+1; Ifim = Ifi-1; Ifip = Ifi+1
	Ila = ImaAll -nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG+1; Jfim = Jfi-1; Jfip = Jfi+1
	Jla = JmaAll -nG; Jlap = Jla+1; Jlam = Jla-1

	u[Ifim:Ila, Jfi:Jla] .= 0.5 * (
		rhou[Ifi:Ilap, Jfi:Jla] + rhou[Ifim:Ila, Jfi:Jla]
	)

	v[Ifi:Ila, Jfim:Jla] .= 0.5 * (
		rhov[Ifi:Jla, Jfi:Jlap] + rhov[Ifi:Jla, Jfim:Jla]
	)

	return nothing
end

# ╔═╡ d7b3b52c-5da3-49a8-b727-375392c327f1
let
	rhou = ones(12,12)
	rhov = ones(12,12)
	u = zeros(12,12)
	v = zeros(12,12)
	mom2vel!(u,v,rhou,rhov)
	u
end

# ╔═╡ 2e300000-a24d-4c5b-b83e-90b1d6eaaa42
function calcFluxConCDS!(fluxConX, fluxConY, Phi, U, V, Dx)

	# using central differencing
	
	nG = 1

	ImaAll, JmaAll = size(fluxConX)
	Ifi = nG+1; Ifim = Ifi-1; Ifip = Ifi+1
	Ila = ImaAll -nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG+1; Jfim = Jfi-1; Jfip = Jfi+1
	Jla = JmaAll -nG; Jlap = Jla+1; Jlam = Jla-1

	PhiC = @view Phi[Ifi:Ila, Jfi:Jla]
	PhiE = @view Phi[Ifip:Ilap, Jfi:Jla]; PhiW = @view Phi[Ifim:Ilam, Jfi:Jla]
	PhiN = @view Phi[Ifi:Ila, Jfip:Jlap]; PhiS = @view Phi[Ifi:Ila, Jfim:Jlam]

	Ue = @view U[Ifi:Ila, Jfi:Jla]; Uw = @view U[Ifim:Ilam, Jfi:Jla]
	Vn = @view V[Ifi:Ila, Jfi:Jla]; Vs = @view V[Ifi:Ila, Jfim:Jlam]

	fluxConX[Ifi:Ila, Jfi:Jla] .= 0.5 *Dx^2 *((PhiE +PhiC) .*Ue .- (PhiC +PhiW) .*Uw)

	fluxConY[Ifi:Ila, Jfi:Jla] .= 0.5 *Dx^2 *((PhiN +PhiC) .*Vn .- (PhiC + PhiS) .*Vs)
	
	return nothing
end

# ╔═╡ c982c464-c55a-49d3-bcea-4690ccf6f28e
let
	Dx = 1 #fluxConX, fluxConY, Phi, U, V,
	Phi = ones(12,12)
	U = ones(12,12)
	V = ones(12,12)
	fluxConX = zeros(12,12)
	fluxConY = zeros(12,12)
	calcFluxConCDS!(fluxConX, fluxConY, Phi, U, V, Dx)
end

# ╔═╡ 49eaf1e1-417e-4978-b596-d82697d993cc
function calcFluxDif!(fluxDifX, fluxDifY, Phi, D, deltaX)

	#use central differencing
	
	nG = 1

	ImaAll, JmaAll = size(fluxDifX)
	Ifi = nG+1; Ifim = Ifi-1; Ifip = Ifi+1
	Ila = ImaAll -nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG+1; Jfim = Jfi-1; Jfip = Jfi+1
	Jla = JmaAll -nG; Jlap = Jla+1; Jlam = Jla-1

	PhiC = @view Phi[Ifi:Ila,Jfi:Jla]
	PhiE = @view Phi[Ifip:Ilap, Jfi:Jla]; PhiW = @view Phi[Ifim:Ilam, Jfi:Jla]
	PhiN = @view Phi[Ifi:Ila, Jfip:Jlap]; PhiS = @view Phi[Ifi:Ila, Jfim:Jlam]

	fluxDifX[Ifi:Ila, Jfi:Jla] .= D *deltaX * (PhiE .- 2*PhiC .+ PhiW)
	fluxDifY[Ifi:Ila, Jfi:Jla] .= D *deltaX * (PhiN .- 2*PhiC .+ PhiS)

	return nothing
end

# ╔═╡ eb9db57e-d6c5-4240-8aa8-0aa22e50ea91
let
	Dx = 1.
	D = 1.
	Phi = ones(12,12)
	fluxDifX = zeros(12,12)
	fluxDifY = zeros(12,12)
	calcFluxDif!(fluxDifX, fluxDifY, Phi, D, Dx)
end

# ╔═╡ 1724d8a9-e494-40ab-96a8-30d1685452c4
function calcFluxConUSCDS!(fluxConX, fluxConY, Phi, U, V, deltaT, deltaX)

	#upwind-shifted CDS scheme 
	
	nG = 1

	ImaAll, JmaAll = size(fluxConX)
	Ifi = nG+1; Ifim = Ifi-1; Ifip = Ifi+1
	Ila = ImaAll -nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG+1; Jfim = Jfi-1; Jfip = Jfi+1
	Jla = JmaAll -nG; Jlap = Jla+1; Jlam = Jla-1

	PhiC = @view Phi[Ifi:Ila, Jfi:Jla]
	PhiE = @view Phi[Ifip:Ilap, Jfi:Jla]; PhiW = @view Phi[Ifim:Ilam, Jfi:Jla]
	PhiN = @view Phi[Ifi:Ila, Jfip:Jlap]; PhiS = @view Phi[Ifi:Ila, Jfim:Jlam]

	Ue = @view U[Ifi:Ila, Jfi:Jla]; Uw = @view U[Ifim:Ilam, Jfi:Jla]
	Vn = @view V[Ifi:Ila, Jfi:Jla]; Vs = @view V[Ifi:Ila, Jfim:Jlam]

	Fe = similar(fluxConX)
	Fw = similar(fluxConX)

	Fe[Ifi:Ila, Jfi:Jla] = (1 + Ue *deltaT/deltaX) *PhiC + (1 - Ue *deltaT/deltaX) *PhiE
	Fw[Ifi:Ila, Jfi:Jla] = (1 + Uw *deltaT/deltaX) *PhiW + (1 - Uw *deltaT/deltaX) *PhiC
	
	fluxConX[Ifi:Ila, Jfi:Jla] .= 0.5 *Dx^2 *( Fe - Fw )

	
	Fn = similar(fluxConX)
	Fs = similar(fluxConX)

	Fn[Ifi:Ila, Jfi:Jla] = (1 + Un *deltaT/deltaX) *PhiC + (1 - Un *deltaT/deltaX) *PhiN
	Fs[Ifi:Ila, Jfi:Jla] = (1 + Us *deltaT/deltaX) *PhiS + (1 - Us *deltaT/deltaX) *PhiC

	fluxConY[Ifi:Ila, Jfi:Jla] .= 0.5 *Dx^2 *( Fn - Fs )
	
	return nothing
end

# ╔═╡ 27f5fa2b-a914-429c-abcd-10acc77bb5e3
let
	Dx = 1.
	Dt = 1.
	D = 1.
	Phi = ones(12,12)
	fluxDifX = zeros(12,12)
	fluxDifY = zeros(12,12)
	calcFluxDif!(fluxDifX, fluxDifY, Phi, D, Dx)
end

# ╔═╡ cf557048-37a6-45a8-940a-7699802c0f29
# revisit scalar transport in 1D 
# add UDS, USCDS, TVD



# ╔═╡ a88e86af-8a0f-4332-8ddb-2a8d09c3d6e2
step = @bind step Slider(0:1:2000, default=200)

# ╔═╡ 1b2e91cb-8ac2-41ad-bad0-b32e32bd1d24


# ╔═╡ 6c803a5c-6ed3-476e-8755-9f7f271eeb56
function step!(phiNew::AbstractVector, phi::AbstractVector, u::AbstractVector, DeltaT::Float64, DeltaX::Float64, D::Float64, nG::Int64)

	ImaAll = length(phiNew)
	Ifim = 1; Ifi = nG+1; 
	Ila = ImaAll -nG; Ilap = Ila+1

	phiNew[Ifi:Ila] = phi[Ifi:Ila] .+ DeltaT * begin
		( -0.5/DeltaX * (phi[Ifi+1:Ilap] .* u[Ifi+1:Ilap] .- phi[Ifim:Ila-1] .* u[Ifim:Ila-1]) ) + ( D/DeltaX^2 * (phi[Ifi+1:Ilap] -2*phi[Ifi:Ila] + phi[Ifim:Ila-1]) )
	end
	
	return nothing
end

# ╔═╡ 09fcf9e4-e6d4-432b-9f5e-0eb0d8474ce7
function applyPeriodicBC!(phi::AbstractVector, nG::Int64)
	
	ImaAll = length(phi)
	Ifim = 1; Ifi = nG+1; 
	Ila = ImaAll -nG; Ilap = Ila+1
	
	phi[Ifim] = phi[Ila]
	phi[Ilap]  = phi[Ifi]
	
	return nothing
end

# ╔═╡ 9773366c-4b37-4405-9768-708b03c93572
function stepUDS!(phiNew::AbstractVector, phi::AbstractVector, u::AbstractVector, DeltaT::Float64, DeltaX::Float64, D::Float64, nG::Int64)

	ImaAll = length(phiNew)
	Ifi = nG+1; Ifim = Ifi-1; Ifip = Ifi+1
	Ila = ImaAll -nG; Ilap = Ila+1; Ilam = Ila-1

	phiC = @view phi[Ifi:Ila]
	phiE = @view phi[Ifip:Ilap]
	phiW = @view phi[Ifim:Ilam]
	
	ux = similar(u); Fe = similar(phi); Fw = similar(phi)
	ux[Ifim:Ila] = 0.5 * ( u[Ifi:Ilap] + u[Ifim:Ila] )

	ue = @view ux[Ifi:Ila]
	uw = @view ux[Ifim:Ilam]
	
	Fe = phiC .*ue .*(ue .>= 0) + phiE .*ue .*(ue .< 0)
	Fw = phiW .*uw .*(uw .>= 0) + phiC .*uw .*(uw .< 0)
	Fc = 1 /DeltaX * (Fe - Fw)
	
	#Fc = 0.5/DeltaX * (phiE .* u[Ifi+1:Ilap] .- phiW .* u[Ifim:Ila-1]) #CDS debug
	
	phiNew[Ifi:Ila] = phi[Ifi:Ila] .+ DeltaT * begin
		-Fc + ( D/DeltaX^2 * (phiE -2*phiC + phiW) )
	end
	
	return nothing
end

# ╔═╡ 970093cc-9c52-4ad7-9d60-1bbfc64f29cc
function stepDDS!(phiNew::AbstractVector, phi::AbstractVector, u::AbstractVector, DeltaT::Float64, DeltaX::Float64, D::Float64, nG::Int64)

	ImaAll = length(phiNew)
	Ifi = nG+1; Ifim = Ifi-1; Ifip = Ifi+1
	Ila = ImaAll -nG; Ilap = Ila+1; Ilam = Ila-1

	phiC = @view phi[Ifi:Ila]
	phiE = @view phi[Ifip:Ilap]
	phiW = @view phi[Ifim:Ilam]
	
	ux = similar(u); Fe = similar(phi); Fw = similar(phi)
	ux[Ifim:Ila] = 0.5 * ( u[Ifi:Ilap] + u[Ifim:Ila] )

	ue = @view ux[Ifi:Ila]
	uw = @view ux[Ifim:Ilam]
	
	Fe = phiC .*ue .*(ue .< 0) + phiE .*ue .*(ue .>= 0)
	Fw = phiW .*uw .*(uw .< 0) + phiC .*uw .*(uw .>= 0)
	Fc = 1 /DeltaX * (Fe - Fw)
	
	#Fc = 0.5/DeltaX * (phiE .* u[Ifi+1:Ilap] .- phiW .* u[Ifim:Ila-1]) #CDS debug
	
	phiNew[Ifi:Ila] = phi[Ifi:Ila] .+ DeltaT * begin
		-Fc + ( D/DeltaX^2 * (phiE -2*phiC + phiW) )
	end
	
	return nothing
end

# ╔═╡ 46552ff0-ad52-4060-aa1e-0daa0a094b66
function stepUSCDS!(phiNew::AbstractVector, phi::AbstractVector, u::AbstractVector, DeltaT::Float64, DeltaX::Float64, D::Float64, nG::Int64)

	ImaAll = length(phiNew)
	Ifi = nG+1; Ifim = Ifi-1; Ifip = Ifi+1
	Ila = ImaAll -nG; Ilap = Ila+1; Ilam = Ila-1

	phiC = @view phi[Ifi:Ila]
	phiE = @view phi[Ifip:Ilap]
	phiW = @view phi[Ifim:Ilam]
	
	ux = similar(u); Fe = similar(phi); Fw = similar(phi)
	ux[Ifim:Ila] = 0.5 * ( u[Ifi:Ilap] + u[Ifim:Ila] )

	ue = @view ux[Ifi:Ila]
	uw = @view ux[Ifim:Ilam]
	
	Fe = (1 .+ ue*DeltaT/DeltaX) .*phiC .+ (1 .- ue*DeltaT/DeltaX) .*phiE
	Fw = (1 .+ uw*DeltaT/DeltaX) .*phiW .+ (1 .- uw*DeltaT/DeltaX) .*phiC
	Fc = 0.5 /DeltaX * (Fe - Fw)
	
	#Fc = 0.5/DeltaX * (phiE .* u[Ifi+1:Ilap] .- phiW .* u[Ifim:Ila-1]) #CDS debug
	
	phiNew[Ifi:Ila] = phi[Ifi:Ila] .+ DeltaT * begin
		-Fc + ( D/DeltaX^2 * (phiE -2*phiC + phiW) )
	end
	
	return nothing
end

# ╔═╡ 2f51bd01-0518-45f0-8143-e6d2ceb8cbca
# 2. Scalar transport in 1D

let
	U = 1.; CFL = 0.25; D = 1e-5
	x = 1.
	
	Ima = 20; nG = 1
	Ifim = 1; Ifi = nG +1; 
	Ila = Ima +nG; Ilap = Ila + 1
	
	u 	= Vector{Float64}(undef, Ilap)
	phi	= similar(u)
	u 	.= U

	println("x = $x, U = $U, CFL = $(CFL), D = $(D)")
	println("Ima = $(Ima)")#, Ifim=$(Ifim), Ifi=$(Ifi), Ila=$(Ila), Ilap=$(Ilap)")

	DeltaX = x / Ima
	DeltaT = CFL *DeltaX /U
	println("DeltaT = $(DeltaT)")
	println("tStep = $(step)")
	
	phi[Ifi:Ila] .= LinRange(DeltaX/2,x - DeltaX/2,Ima) .|> t->sin(2π*t)
	
	applyPeriodicBC!(phi, nG)

	phiNew = similar(phi)

	plot(LinRange(0,x,Ima), phi[Ifi:Ila], legend=nothing, ylims=(-1.2, 1.2), xlims=(0,x), title="time = $(step*DeltaT)s")
	for _ in 1:step
		stepUSCDS!(phiNew, phi, u, DeltaT, DeltaX, D, nG)
		applyPeriodicBC!(phiNew, nG)
		phi = phiNew
		#plot!(LinRange(0,x,Ima), phiNew[Ifi:Ila])
	end
	plot!(LinRange(0,x,Ima), phiNew[Ifi:Ila])
	plot!()
end

# ╔═╡ f36b351a-2287-40c9-8b46-9a8f56c908f2
function stepTVD!(phiNew::AbstractVector, phi::AbstractVector, u::AbstractVector, DeltaT::Float64, DeltaX::Float64, D::Float64, nG::Int64)

	ImaAll = length(phiNew)
	Ifi = nG+1; Ifim = Ifi-1; Ifip = Ifi+1
	Ila = ImaAll -nG; Ilap = Ila+1; Ilam = Ila-1

	phiC = @view phi[Ifi:Ila]
	phiE = @view phi[Ifip:Ilap]
	phiW = @view phi[Ifim:Ilam]

	uC = @view u[Ifi:Ila]
	uE = @view u[Ifip:Ilap]
	uW = @view u[Ifim:Ilam]

	r = similar(u); B = similar(u)
	r[Ifi:Ila] = (phiE -phiC) ./(phiC -phiW)
	BFunc(r) = max(0., min(2*r, 0.75*r +0.25, 4.))#(1.5*r + 0.5)/2 	# could be one of the input
	B[Ifi:Ila] = r[Ifi:Ila] .|> BFunc

	Fc = 0.5 /DeltaX * (
		(phiC .+ B[Ifi:Ila] .*(phiC .-phiW)) .*(uE .+ uC) .- phiW .*(uW .+ uC)
	)
	
	phiNew[Ifi:Ila] = phi[Ifi:Ila] .+ DeltaT * begin
		-Fc + ( D/DeltaX^2 * (phiE -2*phiC + phiW) )
	end
	
	return nothing
end

# ╔═╡ Cell order:
# ╠═7a635be6-ff7c-11ee-2cdd-afdcd8ccc93d
# ╠═718f1af5-c5ea-4b78-b748-eae0c996cf97
# ╠═ad19c163-839b-47af-a528-4df2499a8bd6
# ╠═d7b3b52c-5da3-49a8-b727-375392c327f1
# ╠═2e300000-a24d-4c5b-b83e-90b1d6eaaa42
# ╠═c982c464-c55a-49d3-bcea-4690ccf6f28e
# ╠═49eaf1e1-417e-4978-b596-d82697d993cc
# ╠═eb9db57e-d6c5-4240-8aa8-0aa22e50ea91
# ╠═1724d8a9-e494-40ab-96a8-30d1685452c4
# ╠═27f5fa2b-a914-429c-abcd-10acc77bb5e3
# ╠═cf557048-37a6-45a8-940a-7699802c0f29
# ╟─a88e86af-8a0f-4332-8ddb-2a8d09c3d6e2
# ╠═2f51bd01-0518-45f0-8143-e6d2ceb8cbca
# ╠═1b2e91cb-8ac2-41ad-bad0-b32e32bd1d24
# ╠═6c803a5c-6ed3-476e-8755-9f7f271eeb56
# ╠═09fcf9e4-e6d4-432b-9f5e-0eb0d8474ce7
# ╠═9773366c-4b37-4405-9768-708b03c93572
# ╠═970093cc-9c52-4ad7-9d60-1bbfc64f29cc
# ╠═46552ff0-ad52-4060-aa1e-0daa0a094b66
# ╠═f36b351a-2287-40c9-8b46-9a8f56c908f2

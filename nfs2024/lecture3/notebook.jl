### A Pluto.jl notebook ###
# v0.19.42

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

# ╔═╡ cd8ba894-feb1-11ee-0543-e3e7557ba55c
import Pkg; Pkg.activate("..")

# ╔═╡ 4b5d953b-3b40-4ee7-ab14-fa604cc0a706
using PlutoUI, Plots, LinearAlgebra

# ╔═╡ d2a36a87-b306-45f4-a8dd-86819b388ac9
step = @bind step Slider(0:1:2000, default=200)

# ╔═╡ f3a193b8-0f48-4650-8c3f-d65676150e3d
function step!(phiNew::AbstractVector, phi::AbstractVector, u::AbstractVector, DeltaT::Float64, DeltaX::Float64, D::Float64, nG::Int64)

	ImaAll = length(phiNew)
	Ifim = 1; Ifi = nG+1; 
	Ila = ImaAll -nG; Ilap = Ila+1

	phiNew[Ifi:Ila] = phi[Ifi:Ila] .+ DeltaT * begin
		( -0.5/DeltaX * (phi[Ifi+1:Ilap] .* u[Ifi+1:Ilap] .- phi[Ifim:Ila-1] .* u[Ifim:Ila-1]) ) + ( D/DeltaX^2 * (phi[Ifi+1:Ilap] -2*phi[Ifi:Ila] + phi[Ifim:Ila-1]) )
	end
	return nothing
end

# ╔═╡ a5ed4acd-b91b-4e3f-9bd1-1f19bc076144
function applyPeriodicBC!(phi::AbstractVector, Ifim::Int64, Ifi::Int64, Ila::Int64, Ilap::Int64)
	phi[Ifim] = phi[Ila]
	phi[Ilap]  = phi[Ifi]
	return nothing
end

# ╔═╡ 69fe07df-f85a-4d96-90b6-c6ad44427cab
meshgrid(xs, ys) = [xs[i] for i in 1:length(xs), j in length(ys)], [ys[j] for i in 1:length(xs), j in 1:length(ys)]

# ╔═╡ 69aed578-c1f8-40d6-9461-9c9ab35ba5ec
let
	x = 1:3
	y = 1:5
	(x = x * ones(5)', y = ones(3) * y')
end

# ╔═╡ c95dc54c-b49a-4a37-8ec5-57410cc2e2d8
#quiver(X,Y,quiver=(U*1e-5,V*1e-5))

# ╔═╡ 73638cc6-1a4d-457e-9fd6-9d26d3c29398
# 3. Scalar transport in 2D

begin
	# case setup
	Ima = 100; Jma = 100; 
	delta = 1/Ima; nG = 2

	Ifi = nG +1
	Ila = Ima +nG
	ImaAll = Ima + 2*nG
	Jfi = nG +1
	Jla = Jma +nG;
	JmaAll = Jma + 2*nG

	DS = 1/4
	CFL = 0.1
	tMax = 0.002
	D = 1e-5
	
	# velocity field: Vortex 
	SF = 10.
	x = LinRange( ((-ImaAll/2) +0.5) *delta, ((ImaAll/2) -0.5) *delta, ImaAll )
	y = LinRange( ((-JmaAll/2) +0.5) *delta, ((JmaAll/2) -0.5) *delta, JmaAll )
	Y = ones(length(x)) * y'
	X = x * ones(length(y))'
	R = sqrt.(X.^2 + Y.^2)
	MV = SF./R
	U = MV .* (Y ./R)
	V = MV .* (-X ./R)

	DeltaT = CFL *delta /(maximum(U))
	println("DeltaT = $(DeltaT)")
	
	# initialize scalar field
	Phi = Matrix{Float64}(undef, ImaAll, JmaAll)
	for i in Ifi:Ila
		for j in Jfi:Jla
			Phi[i,j] = sin((i-nG-0.5)/Ima*2.0*pi) * sin((j-nG-0.5)/Jma*2.0*pi)
		end
	end

	Phie = zeros(size(Phi)...)
	Phin = zeros(size(Phi)...)
	Ue   = zeros(size(Phi)...)
	Vn   = zeros(size(Phi)...)
	Fc   = zeros(size(Phi)...)
	Fd   = zeros(size(Phi)...)
	PhiNew = zeros(size(Phi)...)

	nothing
end

# ╔═╡ 4cc27adc-2501-4de8-94cc-76a2326379b1
#Center2Surface!(Ue, Vn, U, V)
#Center2Surface!(Phie, Phin, Phi, Phi)
#ConvFlux!(Fc, Phie, Phin, Ue, Vn, delta, Ifim, Ifi, Ila, Ilap, Jfim, Jfi, Jla, Jlap)
#DiffFlux!(Fd, Phi, delta, D, Ifim, Ifi, Ila, Ilap, Jfim, Jfi, Jla, Jlap)
#step!(PhiNew, Phi, Fc, Fd, DeltaT, Ifi, Ila, Jfi, Jla)

# ╔═╡ 0fce36c6-552e-4e83-a6e6-1cac9e478431
function Center2Surface!(Ue::AbstractMatrix, Vn::AbstractMatrix, Uc::AbstractMatrix, Vc::AbstractMatrix, nG::Int64)

	ImaAll, JmaAll = size(Uc)
	Ifi = nG+1; Ifip = Ifi+1; Ifim = Ifi-1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Ilapp = Ilap+1; Ifimm = Ifim-1
	Jfi = nG+1; Jfip = Jfi+1; Jfim = Jfi-1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1
	Jlapp = Jlap+1; Jfimm = Jfim-1
	
	Ue[Ifim:Ila, Jfi:Jla] = 0.5 *(Uc[Ifi:Ilap, Jfi:Jla] + Uc[Ifim:Ila, Jfi:Jla])
	Vn[Ifi:Ila, Jfim:Jla] = 0.5 *(Vc[Ifi:Ila, Jfi:Jlap] + Vc[Ifi:Ila, Jfim:Jla])

	#Ue[1:end-1, 2:end-1] .= 0.5 * (Uc[2:end, 2:end-1] + Uc[1:end-1, 2:end-1])
	#Vn[2:end-1, 1:end-1] .= 0.5 * (Vc[2:end-1, 2:end] + Vc[2:end-1, 1:end-1])

	#Ue[end, :] .= Ue[1, :]
	#Vn[:, end] .= Vc[:, 1]
	
	return nothing
end

# ╔═╡ 91b5eb92-01da-4d6a-9f30-7615d7dc8673
function applyPeriodicBC!(Phi::AbstractMatrix, nG::Int64)

	ImaAll, JmaAll = size(Phi)
	Ifi = nG+1; Ifip = Ifi+1; Ifim = Ifi-1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG+1; Jfip = Jfi+1; Jfim = Jfi-1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1
	
	Phi[Ifi:Ila, Jfim] .= Phi[Ifi:Ila, Jla]
	Phi[Ifi:Ila, Jlap] .= Phi[Ifi:Ila, Jfi]
	Phi[Ifim, Jfi:Jla] .= Phi[Ila, Jfi:Jla]
	Phi[Ilap, Jfi:Jla] .= Phi[Ifi, Jfi:Jla]
	
	#Phi[Ifi:Ila, Jfim-1] .= Phi[Ifi:Ila, Jfim]
	#Phi[Ifi:Ila, Jlap+1] .= Phi[Ifi:Ila, Jlap]
	#Phi[Ifim-1, Jfi:Jla] .= Phi[Ifim, Jfi:Jla]
	#Phi[Ilap+1, Jfi:Jla] .= Phi[Ilap, Jfi:Jla]
	
	return nothing
end

# ╔═╡ ecddc9e5-9523-48ba-afd3-394c27aaeb47
function ConvFlux!(Fc::AbstractMatrix, Phie::AbstractMatrix, Phin::AbstractMatrix, Ue::AbstractMatrix, Vn::AbstractMatrix, delta::Float64, nG::Int64)

	ImaAll, JmaAll = size(Fc)
	Ifi = nG+1; Ifip = Ifi+1; Ifim = Ifi-1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG+1; Jfip = Jfi+1; Jfim = Jfi-1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1

	Phinn = @view Phin[Ifi:Ila, Jfi:Jla]; Phins = @view Phin[Ifi:Ila, Jfim:Jlam]
	Vnn   = @view   Vn[Ifi:Ila, Jfi:Jla]; Vns   = @view   Vn[Ifi:Ila, Jfim:Jlam]
	
	Phiee = @view Phie[Ifi:Ila, Jfi:Jla]; Phiew = @view Phie[Ifim:Ilam, Jfi:Jla]
	Uee   = @view   Ue[Ifi:Ila, Jfi:Jla]; Uew   = @view   Ue[Ifim:Ilam, Jfi:Jla]
	
	Fc[Ifi:Ila, Jfi:Jla] = 1 /delta * (

		Phinn .*Vnn - Phins .*Vns + Phiee .*Uee - Phiew .*Uew
		
	)
	return nothing
end

# ╔═╡ e9ca7658-82cb-4e9d-9fce-3edc6d1635aa
function DiffFlux!(Fd::AbstractMatrix, Phi::AbstractMatrix, delta::Float64, D::Float64, nG::Int64)

	ImaAll, JmaAll = size(Fc)
	Ifi = nG+1; Ifip = Ifi+1; Ifim = Ifi-1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG+1; Jfip = Jfi+1; Jfim = Jfi-1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1

	PhiC = @view Phi[Ifi:Ila, Jfi:Jla]
	PhiN = @view Phi[Ifip:Ilap, Jfi:Jla]
	PhiS = @view Phi[Ifim:Ilam, Jfi:Jla]
	PhiE = @view Phi[Ifi:Ila, Jfip:Jlap]
	PhiW = @view Phi[Ifi:Ila, Jfim:Jlam]
	
	Fd[Ifi:Ila, Jfi:Jla] = D /delta^2 * (

		-4 *PhiC + PhiN + PhiS + PhiE + PhiW
		
	)
	return nothing
end

# ╔═╡ 19468266-b1d8-457c-ac10-df9ec3a06b3a
function step!(PhiNew::AbstractMatrix, Phi::AbstractMatrix, Fc::AbstractMatrix, Fd::AbstractMatrix, DeltaT::Float64, nG::Int64)

	ImaAll, JmaAll = size(Fc)
	Ifi = nG+1
	Ila = ImaAll - nG
	Jfi = nG+1
	Jla = JmaAll - nG
	
	PhiNew[Ifi:Ila, Jfi:Jla] = Phi[Ifi:Ila, Jfi:Jla] .+ DeltaT * (
		
		- Fc[Ifi:Ila, Jfi:Jla] + Fd[Ifi:Ila, Jfi:Jla]
		
	)
	
	return nothing
end

# ╔═╡ 5069abb7-9f54-4390-83cd-847ac185bee7
# 2. Scalar transport in 1D

let
	U = 1.; CFL = 0.25; D = 1e-5
	x = 1.
	
	Ima = 100; nG = 1
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
	
	applyPeriodicBC!(phi, Ifim, Ifi, Ila, Ilap)

	phiNew = similar(phi)

	plot(LinRange(0,x,Ima), phi[Ifi:Ila], legend=nothing, ylims=(-1.2, 1.2), xlims=(0,x), title="time = $(step*DeltaT)s")
	for _ in 1:step
		step!(phiNew, phi, u, DeltaT, DeltaX, D, nG)
		applyPeriodicBC!(phiNew, Ifim, Ifi, Ila, Ilap)
		phi = phiNew
		#plot!(LinRange(0,x,Ima), phiNew[Ifi:Ila])
	end
	plot!(LinRange(0,x,Ima), phiNew[Ifi:Ila])
	plot!()
end

# ╔═╡ 9be23c76-87f6-423f-a04b-719d5b4b036f
function simulate!(PhiNew::AbstractMatrix, Fc::AbstractMatrix, Fd::AbstractMatrix, Phie::AbstractMatrix, Phin::AbstractMatrix, Ue::AbstractMatrix, Vn::AbstractMatrix, Phi::AbstractMatrix, U::AbstractMatrix, V::AbstractMatrix, delta::Float64, D::Float64, DeltaT::Float64, tMax::Float64, nG::Int64)
	
	applyPeriodicBC!(Phi, nG)
	applyPeriodicBC!(U, nG); applyPeriodicBC!(V, nG)
	
	Center2Surface!(Ue, Vn, U, V, nG)
	applyPeriodicBC!(Ue, nG); applyPeriodicBC!(Vn, nG)
	
	T = 0.
	while T < tMax
		Center2Surface!(Phie, Phin, Phi, Phi, nG)
		#applyPeriodicBC!(Phie, nG); applyPeriodicBC!(Phin, nG)
		ConvFlux!(Fc, Phie, Phin, Ue, Vn, delta, nG)
		DiffFlux!(Fd, Phi, delta, D, nG)
		step!(PhiNew, Phi, Fc, Fd, DeltaT, nG)
		applyPeriodicBC!(PhiNew, nG)
		Phi = PhiNew
		T += DeltaT
	end
	return nothing
end

# ╔═╡ b01cec0b-bc64-4221-85b9-f257e879f229
@time simulate!(PhiNew, Fc, Fd, Phie, Phin, Ue, Vn, Phi, U, V, delta, D, DeltaT, tMax, nG)

# ╔═╡ d1f9e059-4020-4f74-a59f-343e9f30f125
heatmap(PhiNew[Ifi:Ila, Jfi:Jla], title="T = $(tMax)s", colorbar=false, size=(750,700))

# ╔═╡ d163e3c7-62c7-45c2-9bf9-191015f86fc5
i = @bind i Slider(1:1:10000, default=1)

# ╔═╡ e16484be-d8e3-4bc0-a12b-ef8ecdabfb9e
# ╠═╡ disabled = true
#=╠═╡
@time begin
	applyPeriodicBC!(Phi, nG)
	applyPeriodicBC!(U, nG); applyPeriodicBC!(V, nG)
	
	Center2Surface!(Ue, Vn, U, V, nG)
	applyPeriodicBC!(Ue, nG); applyPeriodicBC!(Vn, nG)
	

	for _ in 1:i
		Center2Surface!(Phie, Phin, Phi, Phi, nG)
		applyPeriodicBC!(Phie, nG); applyPeriodicBC!(Phin, nG)
		ConvFlux!(Fc, Phie, Phin, Ue, Vn, delta, nG)
		DiffFlux!(Fd, Phi, delta, D, nG)
		step!(PhiNew, Phi, Fc, Fd, DeltaT, nG)
		applyPeriodicBC!(PhiNew, nG)
		Phi = copy(PhiNew)
	end

	#heatmap(PhiNew[Ifi:Ila, Jfi:Jla], title="$(i*DeltaT)s")
end
  ╠═╡ =#

# ╔═╡ d518331f-eda4-4d28-ba53-e054b3b7546a
function simulate2!(PhiNew::AbstractMatrix, Fc::AbstractMatrix, Fd::AbstractMatrix, Phie::AbstractMatrix, Phin::AbstractMatrix, Ue::AbstractMatrix, Vn::AbstractMatrix, Phi::AbstractMatrix, U::AbstractMatrix, V::AbstractMatrix, delta::Float64, D::Float64, DeltaT::Float64, stepMax::Int64, nG::Int64)
	
	applyPeriodicBC!(Phi, nG)
	applyPeriodicBC!(U, nG); applyPeriodicBC!(V, nG)
	
	Center2Surface!(Ue, Vn, U, V, nG)
	Center2Surface!(Phie, Phin, Phi, Phi, nG)
	
	applyPeriodicBC!(Ue, nG); applyPeriodicBC!(Vn, nG)
	applyPeriodicBC!(Phie, nG); applyPeriodicBC!(Phin, nG)

	res = Vector{typeof(Phi)}(undef, stepMax+1)
	step = 1
	res[1] = Phi
	T = 0.
	while step < stepMax
		ConvFlux!(Fc, Phie, Phin, Ue, Vn, delta, nG)
		DiffFlux!(Fd, Phi, delta, D, nG)
		step!(PhiNew, Phi, Fc, Fd, DeltaT, nG)
		applyPeriodicBC!(PhiNew, nG)
		Phi = PhiNew
		T += DeltaT
		step += 1
		res[step] = copy(PhiNew)
	end
	return res
end

# ╔═╡ 6239fbd6-0d08-42a8-9496-6c732a119253
begin
	#PhiNew2 = simulate2!(PhiNew, Fc, Fd, Phie, Phin, Ue, Vn, Phi, U, V, delta, D, DeltaT, 1000, nG)
end

# ╔═╡ a4a7d30a-0c06-4ba7-a5df-5cf352e76b6b
#i = @bind i Slider(1:1:1000, default=1)

# ╔═╡ 993e4a10-a173-4302-9163-0956429e4497
#surface(PhiNew2[i])

# ╔═╡ Cell order:
# ╠═cd8ba894-feb1-11ee-0543-e3e7557ba55c
# ╠═4b5d953b-3b40-4ee7-ab14-fa604cc0a706
# ╠═d2a36a87-b306-45f4-a8dd-86819b388ac9
# ╠═5069abb7-9f54-4390-83cd-847ac185bee7
# ╠═f3a193b8-0f48-4650-8c3f-d65676150e3d
# ╟─a5ed4acd-b91b-4e3f-9bd1-1f19bc076144
# ╟─69fe07df-f85a-4d96-90b6-c6ad44427cab
# ╠═69aed578-c1f8-40d6-9461-9c9ab35ba5ec
# ╠═c95dc54c-b49a-4a37-8ec5-57410cc2e2d8
# ╠═73638cc6-1a4d-457e-9fd6-9d26d3c29398
# ╟─4cc27adc-2501-4de8-94cc-76a2326379b1
# ╠═0fce36c6-552e-4e83-a6e6-1cac9e478431
# ╠═91b5eb92-01da-4d6a-9f30-7615d7dc8673
# ╠═ecddc9e5-9523-48ba-afd3-394c27aaeb47
# ╠═e9ca7658-82cb-4e9d-9fce-3edc6d1635aa
# ╠═19468266-b1d8-457c-ac10-df9ec3a06b3a
# ╠═9be23c76-87f6-423f-a04b-719d5b4b036f
# ╠═b01cec0b-bc64-4221-85b9-f257e879f229
# ╠═d1f9e059-4020-4f74-a59f-343e9f30f125
# ╠═d163e3c7-62c7-45c2-9bf9-191015f86fc5
# ╟─e16484be-d8e3-4bc0-a12b-ef8ecdabfb9e
# ╟─d518331f-eda4-4d28-ba53-e054b3b7546a
# ╟─6239fbd6-0d08-42a8-9496-6c732a119253
# ╠═a4a7d30a-0c06-4ba7-a5df-5cf352e76b6b
# ╠═993e4a10-a173-4302-9163-0956429e4497

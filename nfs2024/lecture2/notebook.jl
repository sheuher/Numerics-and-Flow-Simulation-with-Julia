### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 4087e45d-48af-4f75-ae6f-ad015527a276
import Pkg; Pkg.activate("..")

# ╔═╡ 8b8b5bea-43c5-4c34-a120-5082e97892e3
using Plots, LinearAlgebra

# ╔═╡ b3e73488-fe36-11ee-1357-253ec5473df0
# 1. Interpolation in 1D

let
	X = 0:10
	Y = @. 3X^3 + 2X^2 + X
	Xmid = 1/2 *( X[1:end-1] + X[2:end] )
	Ymid = 1/2 *( Y[1:end-1] + Y[2:end] )
	scatter(X,Y, label="Y")
	scatter!(Xmid, Ymid, label="Ymid")

	delta = 0.1
	Xf = 0:delta:10
	dX = 0:delta:1
	Yfine = similar(Xf)
	j = 0
 	for i in 1:length(X)-1
		for d in 1:length(dX)
			Yfine[d + (i-1)*10] = dX[d]/delta *Y[i+1] + (delta-dX[d])/delta *Y[i]
		end
	end

	
	#scatter!(Xf, Yfine, alpha=0.1, label="Yfine")
	scatter(Xf, Yfine)
	scatter!(X,Y)
	scatter!(Xmid,Ymid)
end


# ╔═╡ 4fe077ae-177e-474c-9ec8-58e07d1444f5
function Center2Surface(Uc::AbstractMatrix, Vc::AbstractMatrix)
	@assert size(Uc) == size(Vc)
	
	Ue = zeros(size(Uc)...)
	Vn = zeros(size(Vc)...)
	
	Ue[:, 2:end] .= 0.5 *( Uc[:, 2:end] + Uc[:, 1:end-1] )
	Vn[1:end-1, :] .= 0.5 *( Vc[2:end, :] + Vc[1:end-1, :] )

	(Ue=Ue, Vn=Vn)
end

# ╔═╡ a6d12366-cfda-47b3-8a86-891d61173213
function Center2Surface!(Ue::AbstractMatrix, Vn::AbstractMatrix, Uc::AbstractMatrix, Vc::AbstractMatrix)
	@assert size(Ue) == size(Vn) == size(Vc) == size(Uc)
	
	Ue[:, 2:end] .= 0.5 *( Uc[:, 2:end] + Uc[:, 1:end-1] )
	Vn[1:end-1, :] .= 0.5 *( Vc[2:end, :] + Vc[1:end-1, :] )
end

# ╔═╡ 15c46cf9-3fe1-49a2-bdf1-70747d055fe0
# 2. Interpolation in 2D

let
	Uc = ones(10,10)
	Vc = ones(10,10)
	tmp = 1:10 |> Vector
	tmp = tmp 
	quiver(tmp, tmp, quiver=(Uc, Vc) )

	#Center2Surface(Uc,Vc)

	Ue = zeros(size(Uc)...)
	Vn = zeros(size(Vc)...)
	Center2Surface!(Ue, Vn, Uc, Vc)
	Ue, Vn
end

# ╔═╡ 83e279ce-f119-43b7-8522-de30b2aaa654
### remind parsa of NFS code

# ╔═╡ 2307ec0d-a9bb-419d-b7bf-866b520f53e6
function linearInterpolation(phiW::Float64, phiE::Float64, x::Float64, dx::Float64)
	
	1/dx *(x *phiW + (dx -x) *phiE)
end

# ╔═╡ 4ed2801e-26cc-45c7-b1fe-988b95edbd75
function bilinearInterpolation(phiNW::Float64, phiNE::Float64, phiSE::Float64, phiSW::Float64, x::Float64, y::Float64, dx::Float64, dy::Float64)
	
	phiN = linearInterpolation(phiNW, phiNE, x, dx)
	phiS = linearInterpolation(phiSW, phiSE, x, dx)
	phi  = linearInterpolation(phiN, phiS, y, dy)
	
	phi
end

# ╔═╡ 1eebbe44-115e-445d-b6e0-3c4214a3a983
# 3. Bi-linear interpolation

let
	phiNW = 10.
	phiNE = 5.
	phiSW = 3.
	phiSE = 7.

	#heatmap([phiNW phiNE
	#		 phiSW phiSE])

	N = 10  	# resolution
	dx = 1.
	dy = 1.
	x  = LinRange(0, 1, N )
	y  = LinRange(0, 1, N )

	phi = Matrix{Float64}(undef, N,N)

	for i in 1:N
		for j in 1:N
			phi[i,j] = bilinearInterpolation(phiNW, phiNE, phiSE, phiSW, x[i], y[j], dx, dy)
		end
	end
	plotly()
	surface(phi)

	# test function the x, y , must be linear, but how to test that
	
end

# ╔═╡ b4f6bc51-057f-4380-abcd-747253f05147
begin
	(linearInterpolation(0.,1.,0.5,1.) == 0.5) |> println
	(bilinearInterpolation(1.,1.,1.,1.,0.5,0.5,1.,1.) == 1.) |> println
end

# ╔═╡ fb0345ea-68a7-4537-9ebf-55ef63772a97
# 4. Differentiation in 1D equidistant grids

# with array operation
# with periodic BC treatments

let

	N = 100
	nG = 1
	fim = 1
	fi = fim + nG
	la = fim + N
	lap = la + nG

	X = Vector{Float64}(undef, lap)
	X[fi:la] .= LinRange(-10,10,N)
	X[fim] = X[la]; X[lap] = X[fi]
	
	Y = X .|> x->sin(x)*x^3
	
	#LDS
	dYdX_LDS = similar(X)
	dYdX_LDS[fim:la] .= (Y[fi:lap] - Y[fim:la]) ./ (X[fi:lap] - X[fim:la])
	dYdX_LDS[fim] = dYdX_LDS[la]; dYdX_LDS[lap] = dYdX_LDS[fi]

	#CDS
	dYdX_CDS = similar(X)
	dYdX_CDS[fi:la] .= (Y[fi+1:lap] - Y[fim:la-1]) ./ (X[fi+1:lap] - X[fim:la-1])
	dYdX_CDS[fim] = dYdX_CDS[la]; dYdX_CDS[lap] = dYdX_CDS[fi]

	#RDS
	dYdX_RDS = similar(X)
	dYdX_RDS[fi:lap] .= (Y[fi:lap] - Y[fim:la]) ./ (X[fi:lap] - X[fim:la])
	dYdX_RDS[fim] = dYdX_RDS[la]; dYdX_RDS[lap] = dYdX_RDS[fi]

	p1 = plot(label=nothing)
	plot!(p1,X[fi:la],Y[fi:la], label="y")
	plot!(p1,X[fi:la],dYdX_LDS[fi:la], label="y′ lds")
	plot!(p1,X[fi:la],dYdX_CDS[fi:la], label="y′ cds")
	plot!(p1,X[fi:la],dYdX_RDS[fi:la], label="y′ rds")

	# higher orders
	d2YdX_CDS = similar(X)
	d2YdX_CDS[fi:la] .= (dYdX_CDS[fi+1:lap] - dYdX_CDS[fim:la-1]) ./ (X[fi+1:lap] - X[fim:la-1])
	d2YdX_CDS[fim] = d2YdX_CDS[la]; d2YdX_CDS[lap] = d2YdX_CDS[fi]

	d3YdX_CDS = similar(X)
	d3YdX_CDS[fi:la] .= (d2YdX_CDS[fi+1:lap] - d2YdX_CDS[fim:la-1]) ./ (X[fi+1:lap] - X[fim:la-1])
	d3YdX_CDS[fim] = d3YdX_CDS[la]; d3YdX_CDS[lap] = d3YdX_CDS[fi]

	d4YdX_CDS = similar(X)
	d4YdX_CDS[fi:la] .= (d3YdX_CDS[fi+1:lap] - d3YdX_CDS[fim:la-1]) ./ (X[fi+1:lap] - X[fim:la-1])
	d4YdX_CDS[fim] = d4YdX_CDS[la]; d4YdX_CDS[lap] = d4YdX_CDS[fi]

	p2 = plot(label=nothing,ylims=(-2.5e3,5e3))
	plot!(p2,X[fi:la],Y[fi:la], label="y")
	plot!(p2,X[fi:la],dYdX_CDS[fi:la], label="y′")
	plot!(p2,X[fi:la],d2YdX_CDS[fi:la], label="y″")
	plot!(p2,X[fi:la],d3YdX_CDS[fi:la], label="y‴")
	plot!(p2,X[fi:la],d4YdX_CDS[fi:la], label="y⁗")

	plot(p1,p2,layout=(2,1),size=(750,750))
end

# ╔═╡ 434dfd92-dd63-4899-aaae-610adaf7ee7f
# 5. Differentiation in 1D non-equidistant grids

# with for loop
# with periodic BC treatments

let
	dx1 = 0.1
	nG = 1
	
	X = [-10.]; dx = dx1
	while X[end] < 10
		push!(X, X[end] + dx)
		dx *= 0.9999
	end

	N = length(X)
	fim = 1
	fi = fim + nG
	la = fim + N
	lap = la + nG
	pushfirst!(X, X[end]); push!(X, X[fi])

	Y = X .|> x->sin(x)*x^3

	#LDS
	dYdX_LDS = similar(X)
	for i in fim:la
		dYdX_LDS[i] = (Y[i+1] - Y[i]) ./ (X[i+1] - X[i])
	end
	dYdX_LDS[fim] = dYdX_LDS[la]; dYdX_LDS[lap] = dYdX_LDS[fi]

	#CDS
	dYdX_CDS = similar(X)
	for i in fi:la
		dYdX_CDS[i] = (Y[i+1] - Y[i-1]) ./ (X[i+1] - X[i-1])
	end
	dYdX_CDS[fim] = dYdX_CDS[la]; dYdX_CDS[lap] = dYdX_CDS[fi]

	#RDS
	dYdX_RDS = similar(X)
	for i in fi:lap
		dYdX_RDS[i] = (Y[i] - Y[i-1]) ./ (X[i] - X[i-1])
	end
	dYdX_RDS[fim] = dYdX_RDS[la]; dYdX_RDS[lap] = dYdX_RDS[fi]
	
	p1 = plot(label=nothing)
	plot!(p1,X[fi:la],Y[fi:la], label="y")
	plot!(p1,X[fi:la],dYdX_LDS[fi:la], label="y′ lds")
	plot!(p1,X[fi:la],dYdX_CDS[fi:la], label="y′ cds")
	plot!(p1,X[fi:la],dYdX_RDS[fi:la], label="y′ rds")
	
	# higher orders
	d2YdX_CDS = similar(X)
	for i in fi:la
		d2YdX_CDS[i] = (dYdX_CDS[i+1] - dYdX_CDS[i-1]) ./ (X[i+1] - X[i-1])
	end
	d2YdX_CDS[fim] = d2YdX_CDS[la]; d2YdX_CDS[lap] = d2YdX_CDS[fi]

	d3YdX_CDS = similar(X)
	for i in fi:la
		d3YdX_CDS[i] = (d2YdX_CDS[i+1] - d2YdX_CDS[i-1]) ./ (X[i+1] - X[i-1])
	end
	d3YdX_CDS[fim] = d3YdX_CDS[la]; d3YdX_CDS[lap] = d3YdX_CDS[fi]

	d4YdX_CDS = similar(X)
	for i in fi:la
		d4YdX_CDS[i] = (d3YdX_CDS[i+1] - d3YdX_CDS[i-1]) ./ (X[i+1] - X[i-1])
	end
	d4YdX_CDS[fim] = d4YdX_CDS[la]; d4YdX_CDS[lap] = d4YdX_CDS[fi]

	p2 = plot(label=nothing,ylims=(-2.5e3,5e3))
	plot!(p2,X[fi:la],Y[fi:la], label="y")
	plot!(p2,X[fi:la],dYdX_CDS[fi:la], label="y′")
	plot!(p2,X[fi:la],d2YdX_CDS[fi:la], label="y″")
	plot!(p2,X[fi:la],d3YdX_CDS[fi:la], label="y‴")
	plot!(p2,X[fi:la],d4YdX_CDS[fi:la], label="y⁗")

	plot(p1,p2,layout=(2,1),size=(750,750))
end

# ╔═╡ 37ebbfea-5e17-49d1-bcdc-37f21c8e2a3f
# 6. Integration

let
	N = 100
	nG = 1
	fim = 1
	fi = fim + nG
	la = fim + N
	lap = la + nG

	X = Vector{Float64}(undef, lap)
	X[fi:la] .= LinRange(-10,10,N)
	X[fim] = X[la]; X[lap] = X[fi]
	filter!(x->3<=x<=9, X)
	Y = X .|> x->sin(x)*x^3

	I = (X[2]-X[1]) *sum(Y)
end

# ╔═╡ 6f647f2b-d66b-4662-a098-a3d21df46e8f
let
	dx1 = 0.1
	nG = 1
	
	X = [-10.]; dx = dx1
	while X[end] < 10
		push!(X, X[end] + dx)
		dx *= 0.9999
	end

	N = length(X)
	fim = 1
	fi = fim + nG
	la = fim + N
	lap = la + nG
	pushfirst!(X, X[end]); push!(X, X[fi])
	filter!(x->3<=x<=9, X)
	Y = X .|> x->sin(x)*x^3

	Dx = 0.5 *(X[3:end] - X[1:end-2])
	push!(Dx, Dx[end])
	pushfirst!(Dx, Dx[1])
	
	I = sum(Y .* Dx)
end

# ╔═╡ ea3528f3-f68b-45cb-9c01-67dea2b22234
# 7. Integration

let
	#Lx = 1
	#Ly = 1
	#A = ∱ dA
	#dA = w(y) *dy
	#w(y) = 1 - y

	y = LinRange(0,1,101)
	dy = y[2] - y[1]
	w = 1 .- y
	A = dy *sum(w)
	println("Area = $(A)")

	yCG = 1/A * dy * sum(y .*w)
	println("yCG = $(yCG)")

	Ixx = 1/A * dy * sum(y .^2 .*w)
	println("Ixx = $(Ixx)")
end

# ╔═╡ Cell order:
# ╠═4087e45d-48af-4f75-ae6f-ad015527a276
# ╠═8b8b5bea-43c5-4c34-a120-5082e97892e3
# ╠═b3e73488-fe36-11ee-1357-253ec5473df0
# ╠═15c46cf9-3fe1-49a2-bdf1-70747d055fe0
# ╠═4fe077ae-177e-474c-9ec8-58e07d1444f5
# ╠═a6d12366-cfda-47b3-8a86-891d61173213
# ╠═1eebbe44-115e-445d-b6e0-3c4214a3a983
# ╠═83e279ce-f119-43b7-8522-de30b2aaa654
# ╠═2307ec0d-a9bb-419d-b7bf-866b520f53e6
# ╠═4ed2801e-26cc-45c7-b1fe-988b95edbd75
# ╠═b4f6bc51-057f-4380-abcd-747253f05147
# ╠═fb0345ea-68a7-4537-9ebf-55ef63772a97
# ╠═434dfd92-dd63-4899-aaae-610adaf7ee7f
# ╠═37ebbfea-5e17-49d1-bcdc-37f21c8e2a3f
# ╠═6f647f2b-d66b-4662-a098-a3d21df46e8f
# ╠═ea3528f3-f68b-45cb-9c01-67dea2b22234

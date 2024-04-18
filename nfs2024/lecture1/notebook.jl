### A Pluto.jl notebook ###
# v0.19.38

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

# ╔═╡ 679afa60-fcc9-11ee-2f33-ddb09de46901
import Pkg; Pkg.activate("..")

# ╔═╡ 5965d931-de1b-416b-8419-9aba97b51868
using PlutoUI, Plots, LinearAlgebra

# ╔═╡ 46600b97-b0d8-4981-a0cf-36c5989e0ba0
begin 
	function dotProduct(x::Vector{Float64})
		x' * x
	end

	function piMonteCarlo(n::Int64)
		n_circle = 0
		for i in 1:n
			if dotProduct(rand(2)) <= 1
				n_circle += 1
			end
		end
		4 * n_circle /n 
	end

	ns = [1:250:100000]
	res = piMonteCarlo.(ns...)
end

# ╔═╡ d9bfbe80-0e1d-46ce-8bdc-eb55d3af9ca1
let
	p = plot(size(250, 250), label=nothing, ylim=(3.,3.25))
	plot!(p, ns, res, label="monte carlo results")
	hline!(p, [π], label="π")
	title!(p, "11. Monte Carlo Method")

	ϵ = abs.( res .- π )
	
	p2 = plot(size(250, 250), label=nothing, ylim=(3.,3.25))
	plot!(p2, ns, ϵ, label="abs error", ylim=(0,0.1))

	plot(p, p2, layout=(2,1))
end

# ╔═╡ 05ae6fed-3edf-4cce-aef2-51e341ea297f
begin 
	x = [1:300]
	y = rand(300)
	z = Vector{Float64}(undef, 298)
	test = []
	#for i in x[2:end-1]
#		z[i-1] = 0.25  *y[i-1] + 0.5 *y[i] + 0.25 *y[i+1]
	#end
	#z = kron(y, diag([0.25,0.5,0.25]))
	z = 0.25 *y[3:end] + 0.5 *y[2:end-1] + 0.25 *y[1:end-2]
end

# ╔═╡ 69fc64b5-0c5f-4edc-a3bd-9a9a3f4eb3b0
plot(
	plot(y, title="10. Filtering", label=nothing),
	plot(z, label = nothing),
	layout=(2,1),
)

# ╔═╡ 3cbd92e1-9039-4c5d-ba7b-47b62c4f7b60
#9. Matrix and Array Operations

let
	x = LinRange(-10,10,200)
	y = sin.(1 ./x)
	plot(x,y)
end

# ╔═╡ be497665-e424-45e7-9687-d16f10699be5
let
	M = [
		1 2 3 5 8.5 2
		-1/3 -2/3 -1 3 2 -5.1
		7 8 9 2.01 -1 3
	]
	N = M[:,4:end]
end

# ╔═╡ f41b62c7-3248-4a2a-a6be-6d8dca267ca8
#8.2 Poisson equation

begin  
	T = zeros(11,11)
	T[1,  	2:end-1] .= 500.
	T[end, 	2:end-1] .= 400.
	T[2:end-1, 	1]   .= 450.
	T[2:end-1,  end] .= 200.
end

# ╔═╡ 114b2f8d-0a22-474a-a31a-490f190af726
step = @bind step Slider(0:10:1000, default=100)

# ╔═╡ f7a4179f-68f3-4ab8-82d2-47902a8b41d1
begin
	T_new = Matrix{Float64}(undef, 11,11)
	T_new = copy(T)
	for _ in 1:step
		T_new[2:end-1, 2:end-1] .= 1/4 *(T_new[3:end, 2:end-1] + T_new[1:end-2, 2:end-1] + T_new[2:end-1, 3:end] + T_new[2:end-1, 1:end-2])
	end
	heatmap(T_new, title="heat equation with step = $(step)")
end

# ╔═╡ 0d38312d-49d4-4fc6-a38e-3c71379ad6cf
# 7. Matrix and Array Operations

let
	A = [1 2 3 
		 4 5 6 
		 7 8 9]
	B = [2/3-1/5 -1/3 1/3-1/5
		 -1/30 2/3 2/3-3/10
		 -1/3 -1/3 1/3]
	(A=A, B=B, mul=A*B, elmMul=A .*B, invB=inv(B), detA=det(A), detB=det(B))
	
end

# ╔═╡ 380c415f-afe4-4d9b-a711-402fc8eddc6e
# 6. Programming 

let

[println("Hello World!") for _ in 1:10]

x = 1:10 ; y = x.^2 ; plot(x,y)

function myFactorial(n::Int64)
	res = n
	while n != 1
		n -= 1
		res *= n
	end
	res
end

println(myFactorial(4) == factorial(4))

function VSphere(r::Number)
	4/3 *pi * r^3
end

println(VSphere(3.5))
end

# ╔═╡ b164b2ea-3237-45c9-a544-cf2cb271c46d
# 3 - 4. Introduction

let
	p = plot(-5:0.1:5,x->exp(x) - exp(-x))
	p2 = surface( -pi:0.1:pi, -pi:0.1:pi, (x,y)->sin(3x) +cos(3x) )
	
	first = 1.23; second = 2.34; println("first = $(first), second = $(second)")
	tmp = second; first = second; second = first; println("first = $(first), second = $(second)\n")

	firstArray = []
	[push!(firstArray, 10+i) for i in 1:6]
	println("firstArray = $(firstArray)")
	secondArray = [21, 22, 23, 24, 25, 26]
	println("secondArray = $(secondArray)")
	println("size = $(size(secondArray))")
	y = secondArray; x = firstArray
	
	p3 = plot(x, y)
	plot(p, p2, p3, layout=(3,1), size=(600,900))
end

# ╔═╡ 6599167d-cd78-46d5-bf58-f9221d80b066
# 2. Introduction

let
	x = 0.4
	y = [0.1, 0.2, 0.3]

	@. log(x^2 /(1 + x^2 + y))
end

# ╔═╡ Cell order:
# ╠═679afa60-fcc9-11ee-2f33-ddb09de46901
# ╠═5965d931-de1b-416b-8419-9aba97b51868
# ╟─46600b97-b0d8-4981-a0cf-36c5989e0ba0
# ╟─d9bfbe80-0e1d-46ce-8bdc-eb55d3af9ca1
# ╠═05ae6fed-3edf-4cce-aef2-51e341ea297f
# ╟─69fc64b5-0c5f-4edc-a3bd-9a9a3f4eb3b0
# ╠═3cbd92e1-9039-4c5d-ba7b-47b62c4f7b60
# ╠═be497665-e424-45e7-9687-d16f10699be5
# ╠═f41b62c7-3248-4a2a-a6be-6d8dca267ca8
# ╟─114b2f8d-0a22-474a-a31a-490f190af726
# ╠═f7a4179f-68f3-4ab8-82d2-47902a8b41d1
# ╠═0d38312d-49d4-4fc6-a38e-3c71379ad6cf
# ╠═380c415f-afe4-4d9b-a711-402fc8eddc6e
# ╠═b164b2ea-3237-45c9-a544-cf2cb271c46d
# ╠═6599167d-cd78-46d5-bf58-f9221d80b066

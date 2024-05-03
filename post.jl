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

# ╔═╡ 8713debe-0289-11ef-13b3-af514c05498e
import Pkg; Pkg.activate("$(@__DIR__)/..")

# ╔═╡ f5a05494-0544-438c-b7a5-ea882cf0bd30
using PlutoUI, Plots, DelimitedFiles

# ╔═╡ 3fad8888-18f6-42c9-bbfd-3863767c5bdb
include("case3.jl")

# ╔═╡ a1508d3f-d638-4afb-a087-703c451ab9b5
dummy = similar(Phi)

# ╔═╡ 855040bb-9d5c-48a7-85db-849f74c0f7d4
#quiver(X,Y,quiver=(U*1e-7,V*1e-7))

# ╔═╡ f0ab5d2a-e418-4442-a64f-d8dd31968afe
step = @bind step Slider(0:saveNT:nTMax, show_value=true)

# ╔═╡ 72cb1adc-8899-44a1-85ce-b99dc09f9a3b
function datPath(varName, nT)
	savePath = string( rpad(varName, 6, '_'), string(nT, pad=6), ".dat")
	savePath = joinpath(@__DIR__, "dat/$(savePath)")
	if ispath(savePath); return savePath; end
end

# ╔═╡ c0c26222-0588-449b-94a4-ffa8fa85919c
function readDat!(dummy, dat::String)
	open(dat, "r") do file
		dummy .= readdlm(dat)
	end
end

# ╔═╡ 111d6c87-3a39-400b-b960-b375512492c6
begin
	path = datPath("Phi", step)
	println(path)
	readDat!(dummy, path)
	heatmap(dummy[Ifi:Ila, Jfi:Jla], size=(700,700), clims=(-1.2,1.2), c=:Greys		
, cbar=false)
	#savefig(plot!(size=(1500,1500), dpi=750), "taylorGreenVortex.png")
end

# ╔═╡ ce31b927-3bc2-4142-8810-4c52ac401b60
# ╠═╡ disabled = true
#=╠═╡
let
	dummy2 = similar(Phi)
	anim = @animate for i in 0:saveNT:nTMax
		path = datPath("Phi", i)
		readDat!(dummy2, path)
		heatmap(dummy2[Ifi:Ila, Jfi:Jla], size=(700,575), clims=(-1.2,1.2))
	end
	gif(anim)
end
  ╠═╡ =#

# ╔═╡ c2d9a342-ac69-455c-a552-e4890875d3f8
# ╠═╡ disabled = true
#=╠═╡
let
	dummy2 = similar(Phi)
	anim = @animate for i in 0:saveNT:nTMax
		path = datPath("Phi", i)
		readDat!(dummy2, path)
		surface(dummy2[Ifi:Ila, Jfi:Jla], size=(700,575), zlims=(-1.2,1.2), clims=(-1.2,1.2))
	end
	gif(anim)
end
  ╠═╡ =#

# ╔═╡ 64ab2b46-10cb-4a7a-ba6a-b279a3ca3cbe
f=x-> 2*tanh(10x) + 20

# ╔═╡ 25071394-8ae5-4d3a-8a21-e8562a2e0692
f.(-1:0.01:1) |> plot

# ╔═╡ Cell order:
# ╠═8713debe-0289-11ef-13b3-af514c05498e
# ╠═f5a05494-0544-438c-b7a5-ea882cf0bd30
# ╠═3fad8888-18f6-42c9-bbfd-3863767c5bdb
# ╠═a1508d3f-d638-4afb-a087-703c451ab9b5
# ╠═855040bb-9d5c-48a7-85db-849f74c0f7d4
# ╟─f0ab5d2a-e418-4442-a64f-d8dd31968afe
# ╠═111d6c87-3a39-400b-b960-b375512492c6
# ╠═ce31b927-3bc2-4142-8810-4c52ac401b60
# ╠═c2d9a342-ac69-455c-a552-e4890875d3f8
# ╟─72cb1adc-8899-44a1-85ce-b99dc09f9a3b
# ╟─c0c26222-0588-449b-94a4-ffa8fa85919c
# ╠═64ab2b46-10cb-4a7a-ba6a-b279a3ca3cbe
# ╠═25071394-8ae5-4d3a-8a21-e8562a2e0692

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
include("case.jl")

# ╔═╡ a1508d3f-d638-4afb-a087-703c451ab9b5
dummy = similar(p)

# ╔═╡ f0ab5d2a-e418-4442-a64f-d8dd31968afe
step = @bind step Slider(0:saveNT:nTMax, show_value=true)

# ╔═╡ f9addb56-623a-4d52-9d11-319e5c80414e
dummy

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
	path = datPath("rhou", step)
	println(path)
	readDat!(dummy, path)
	heatmap(dummy)
end

# ╔═╡ 34211587-e366-4c5b-8223-c128b0903862
begin
dummyrhov = similar(dummy)
dummyp = similar(dummy)
	1
end

# ╔═╡ c21d103b-3f3d-457b-a501-13a10d200153
begin
	x = LinRange(0,1+delta,ImaAll)
	X = x * ones(length(x))'
	Y = ones(length(x)) * x'
	1
end

# ╔═╡ 2781a5a8-e96c-4a0a-955d-8f1ead79a24e
begin 
	pathrhov = datPath("rhov", step)
	println(pathrhov)
	readDat!(dummyrhov, pathrhov)
	pathp = datPath("p", step)
	readDat!(dummyp, pathp)
	quiver(X, Y, quiver=(dummy*1e-4, dummyrhov*1e-4))
end

# ╔═╡ 0d9248ab-ab45-431c-ad0d-7fef21a956cd
heatmap(dummyp)

# ╔═╡ Cell order:
# ╠═8713debe-0289-11ef-13b3-af514c05498e
# ╠═f5a05494-0544-438c-b7a5-ea882cf0bd30
# ╠═3fad8888-18f6-42c9-bbfd-3863767c5bdb
# ╠═a1508d3f-d638-4afb-a087-703c451ab9b5
# ╟─f0ab5d2a-e418-4442-a64f-d8dd31968afe
# ╠═111d6c87-3a39-400b-b960-b375512492c6
# ╠═f9addb56-623a-4d52-9d11-319e5c80414e
# ╠═72cb1adc-8899-44a1-85ce-b99dc09f9a3b
# ╟─c0c26222-0588-449b-94a4-ffa8fa85919c
# ╠═34211587-e366-4c5b-8223-c128b0903862
# ╠═c21d103b-3f3d-457b-a501-13a10d200153
# ╠═2781a5a8-e96c-4a0a-955d-8f1ead79a24e
# ╠═0d9248ab-ab45-431c-ad0d-7fef21a956cd

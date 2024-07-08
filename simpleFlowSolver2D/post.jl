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
dummy = similar(Phi)

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
	heatmap(dummy[Ifi:Ila, Jfi:Jla])
end

# ╔═╡ Cell order:
# ╠═8713debe-0289-11ef-13b3-af514c05498e
# ╠═f5a05494-0544-438c-b7a5-ea882cf0bd30
# ╠═3fad8888-18f6-42c9-bbfd-3863767c5bdb
# ╠═a1508d3f-d638-4afb-a087-703c451ab9b5
# ╟─f0ab5d2a-e418-4442-a64f-d8dd31968afe
# ╠═111d6c87-3a39-400b-b960-b375512492c6
# ╟─72cb1adc-8899-44a1-85ce-b99dc09f9a3b
# ╟─c0c26222-0588-449b-94a4-ffa8fa85919c

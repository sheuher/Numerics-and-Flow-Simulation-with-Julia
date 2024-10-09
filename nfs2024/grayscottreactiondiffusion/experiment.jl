### A Pluto.jl notebook ###
# v0.19.45

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

# ╔═╡ 1c929afa-8519-4038-8240-0f42c984e64f
begin
	import Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
	using Plots
	using PlutoUI
end

# ╔═╡ 5291c9b2-6396-11ef-06bf-1916151ab9b5
begin
	macro qx(ix, iy, ik) esc(:( (C[$ix+1, $iy, $ik] - C[$ix, $iy, $ik])  *D[$ik] /dx^2 )) end
	macro qy(ix, iy, ik) esc(:( (C[$ix, $iy+1, $ik] - C[$ix, $iy, $ik])  *D[$ik] /dy^2 )) end
	
	macro react(ix, iy, ik) esc(:( 
	                                (-C[$ix, $iy, 1]* C[$ix, $iy, 2] ^2) * ($ik == 1) +
	                                ( C[$ix, $iy, 1]* C[$ix, $iy, 2] ^2) * ($ik == 2)
	                        ))
	end
	
	macro kll(ix, iy, ik) esc(:( 
	                                (F *(1 -C[$ix, $iy, 1])) * ($ik == 1) +
	                                (-(F + k) *C[$ix, $iy, 2]) * ($ik == 2)
	                        ))
	end
	
end

# ╔═╡ 57538fce-15fe-4413-b8c0-b2bd5e85192e
function step!(C2, C, DeltaX, DeltaT)
    nG = 1
    ImaAll, JmaAll, nPhi = size(C)
	Ifi = nG+1;     Ila = ImaAll - nG
	Jfi = nG+1;     Jla = JmaAll - nG

    dx = DeltaX
    dt = DeltaT

    D = [2.0E-05, 1.0E-05]
    F = 0.029
    k = 0.057

    for k in 1:nPhi
        for j in 1:1
            for i in Ifi:Ila
                @inbounds C2[i, j, k] = C[i, j, k] + dt * (
                   #@react(i, j ,k) + @kll(i, j, k)
					C[i,j,1]*C[i,j,k]^2 * -1*(k==1) * 1*(k==2) 
					+ F*(1 - C[i,j,1]) * (k==1)
					- (F + k) * C[i,k,2] * (k==2)
                    + @qx(i, j, k) - @qx(i-1, j, k)
                )
            end
        end
    end
    C2
end

# ╔═╡ b8a18bac-fe4f-4b19-bf75-bc53684e1433
let
		scatter(Phi[:,:,1])
		scatter!(Phi[:,:,2])
end

# ╔═╡ abd106d3-4a9e-4d59-92b2-d1e49a92c6f3
begin
	saveNT = 100
	nTMax = 2000
    D = [2.0E-05, 1.0E-05]
	DeltaT = 0.8 * (2.5/Ima)^2 / 2*(maximum(D)) 
end

# ╔═╡ a5b6c7c8-6a85-4d4c-866f-f8b222dbb19a
	step = @bind step Slider(0:1:nTMax, show_value=true)

# ╔═╡ cecf61f6-116b-4130-bd3c-5eba54815d24
begin
	
Ima = 1000; 
delta = 2.5/Ima; 
nG = 1;
Ifi = nG +1; Ifip=Ifi+1; Ifim=Ifi-1
Ila = Ima +nG; Ilap=Ila+1; Ilam=Ila-1
ImaAll = Ima + 2*nG
nPhi = 2
JmaAll = 1

x = LinRange( ((-ImaAll/2) +0.5) *delta, ((ImaAll/2) -0.5) *delta, ImaAll )
	
Phi  = Array{Float64}(undef, ImaAll, JmaAll, nPhi)
fill!(Phi, 0.0)
Phi[:,:,1] = @. 1 - exp(-80*((x+0.05)^2) )
Phi[:,:,2] = @. exp(-80*((x-0.05)^2))
PhiNew = similar(Phi)
fill!(PhiNew, 0.0)

	for nt = 1:step
		Phi[Ifim, :, :] .= Phi[Ila, :, :]
		Phi[Ilap, :, :] .= Phi[Ifi, :, :]
		global PhiNew = step!(PhiNew, Phi, delta, DeltaT)
		global Phi, PhiNew = PhiNew, Phi
	end
end

# ╔═╡ Cell order:
# ╠═1c929afa-8519-4038-8240-0f42c984e64f
# ╠═5291c9b2-6396-11ef-06bf-1916151ab9b5
# ╠═57538fce-15fe-4413-b8c0-b2bd5e85192e
# ╠═cecf61f6-116b-4130-bd3c-5eba54815d24
# ╠═abd106d3-4a9e-4d59-92b2-d1e49a92c6f3
# ╠═a5b6c7c8-6a85-4d4c-866f-f8b222dbb19a
# ╠═b8a18bac-fe4f-4b19-bf75-bc53684e1433

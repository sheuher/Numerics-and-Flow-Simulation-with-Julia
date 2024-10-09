import Pkg; Pkg.activate("..")
using LinearAlgebra, Plots

include("mom2vel.jl")
include("calcFluxConCDS.jl")
include("calcFluxDif.jl") 
#include("divMom2D.jl") 
#include("poissonSolver2D.jl")
#include("corrMom2D.jl") 
include("case3.jl")
include("applyBC.jl")
include("save.jl")
include("log.jl")
# include("test.jl")

function step1!(PhiNew, Phi, FcX, FcY, FdX, FdY, Ue, Vn, DeltaX, DeltaT)
    nG = 1
    ImaAll, JmaAll = size(Phi)
	Ifi = nG+1;     Ila = ImaAll - nG
	Jfi = nG+1;     Jla = JmaAll - nG

    #applyBC!(Ue); applyBC!(Vn); applyBC!(Phi)

    calcFluxConCDS!(FcX, FcY, Phi, Ue, Vn, DeltaX)

    calcFluxDif!(FdX, FdY, Phi, D, DeltaX)
    
    PhiNew[Ifi:Ila, Jfi:Jla] .= Phi[Ifi:Ila, Jfi:Jla] .+ DeltaT/DeltaX^2 *( 
        .- FcX[Ifi:Ila, Jfi:Jla] .- FcY[Ifi:Ila, Jfi:Jla] .+ FdX[Ifi:Ila, Jfi:Jla] .+ FdY[Ifi:Ila, Jfi:Jla] 
        )
    
    applyBC!(PhiNew)    

    nothing
end

macro qx(ix, iy) esc(:( (C[$ix+1, $iy] - C[$ix, $iy])  *D /dx^2 )) end
macro qy(ix, iy) esc(:( (C[$ix, $iy+1] - C[$ix, $iy])  *D /dy^2 )) end

macro fcx(ix, iy) esc(:( 0.5 /dx * (C[$ix+1, $iy] 
    + C[$ix, $iy]) *Uf[$ix, $iy] )) end
macro fcy(ix, iy) esc(:( 0.5 /dy * (C[$ix, $iy+1]
    + C[$ix, $iy])  *Vf[$ix, $iy] )) end

function step2!(C2, C, FcX, FcY, FdX, FdY, Uf, Vf, DeltaX, DeltaT)
    nG = 1
    ImaAll, JmaAll = size(C)
	Ifi = nG+1;     Ila = ImaAll - nG
	Jfi = nG+1;     Jla = JmaAll - nG

    D = 1e-5
    dx = dy = DeltaX
    dt = DeltaT

    for j in Jfi:Jla
        for i in Ifi:Ila
            @inbounds C2[i, j] = C[i, j] + dt * ( 
                - @fcx(i, j) + @fcx(i-1, j) - @fcy(i, j) + @fcy(i, j-1)
                + @qx(i, j) - @qx(i-1, j) + @qy(i, j) - @qy(i, j-1)
             )
        end
    end

    applyBC!(C2)    

    nothing
end

# t = 0
# I = 0
# mom2vel!(Ue, Vn, U, V)
# save.([Phi,Ue,Vn], ["Phi", "Ue", "Vn"], I)
# log(Ima,Jma,nG,CFL,delta,DeltaT,tMax,nTMax,D,saveNT)
# while (t<tMax && I<=nTMax)
#     step!(PhiNew, Phi, FcX, FcY, FdX, FdY, Ue, Vn, delta, DeltaT)
#     global Phi = PhiNew
#     global t+=DeltaT
#     global I+=1
#     log(I,DeltaT,t)
#     if I % saveNT == 0
#         save.([Phi,Ue,Vn], ["Phi", "Ue", "Vn"], I)
#     end
# end
# log()

function simulate!(PhiNew, Phi, FcX, FcY, FdX, FdY, Ue, Vn, U, V, DeltaX, DeltaT, Ima,Jma,nG,CFL,tMax,nTMax,D,saveNT)
    
    t = 0
    I = 0
    mom2vel!(Ue, Vn, U, V)
    applyBC!(Phi)
    applyBC!(PhiNew)
    applyBC!(Ue)
    applyBC!(Vn)
    save.([Phi,Ue,Vn], ["Phi", "Ue", "Vn"], I)
    log(Ima,Jma,nG,CFL,delta,DeltaT,tMax,nTMax,D,saveNT)
    while (t<tMax && I<=nTMax)
        step2!(PhiNew, Phi, FcX, FcY, FdX, FdY, Ue, Vn, DeltaX, DeltaT)

        Phi, PhiNew = PhiNew, Phi
        t += DeltaT
        I += 1

        log(I,DeltaT,t)
        ( I % saveNT == 0 ) &&
            #save.([Phi,Ue,Vn], ["Phi", "Ue", "Vn"], I)
            save.([Array(Phi)], ["Phi"], I) 
    end
    log()

end

simulate!(PhiNew, Phi, FcX, FcY, FdX, FdY, Ue, Vn, U, V, delta, DeltaT, Ima,Jma,nG,CFL,tMax,nTMax,D,saveNT)
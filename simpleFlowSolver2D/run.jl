import Pkg; Pkg.activate("..")
using LinearAlgebra, Plots
include("mom2vel.jl")
include("calcFluxConCDS.jl")
include("calcFluxDif.jl") 
include("divMom2D.jl") 
include("poissonSolver2D.jl")
include("corrMom2D.jl") 
include("case.jl")
include("applyBC.jl")
include("save.jl")
include("log.jl")

function step!(PhiNew, Phi, FcX, FcY, FdX, FdY, Ue, Vn, U, V, DeltaX, DeltaT)
    nG = 1
    ImaAll, JmaAll = size(Phi)
	Ifi = nG+1
	Ila = ImaAll - nG
	Jfi = nG+1
	Jla = JmaAll - nG
    
    calcFluxConCDS!(FcX, FcY, Phi, Ue, Vn, DeltaX)
    calcFluxDif!(FdX, FdY, Phi, D, DeltaX)

    PhiNew[Ifi:Ila, Jfi:Jla] .= Phi[Ifi:Ila, Jfi:Jla] .+ DeltaT/DeltaX^3 *( -FcX[Ifi:Ila, Jfi:Jla] .- FcY[Ifi:Ila, Jfi:Jla] .+ FdX[Ifi:Ila, Jfi:Jla] .+ FdY[Ifi:Ila, Jfi:Jla] ) 

    applyBC!(PhiNew)

    #divMom2D!(divPred, rhouP, DeltaX)
    #poissonSolver2D!(p, divPred, dt, delta, nItMax, epsMax)
    #corrMom2D!(rhou, rhouP, p, dt, delta)
    mom2vel!(Ue, Vn, U, V)
end

# t = 0
# I = 0
# res = []
# mom2vel!(Ue, Vn, U, V)
# save.([Phi,Ue,Vn], ["Phi", "Ue", "Vn"], I)
# log(Ima,Jma,nG,CFL,delta,DeltaT,tMax,nTMax,D,saveNT)
# while (t<tMax && I<nTMax)
#     log(I,DeltaT,t)
#     step!(PhiNew, Phi, FcX, FcY, FdX, FdY, Ue, Vn, delta, DeltaT)
#     global Phi = PhiNew
#     global t+=DeltaT
#     global I+=1
#     if I % saveNT == 0
#         save.([Phi,Ue,Vn], ["Phi", "Ue", "Vn"], I)
#     end
# end
# log()

function simulate!(PhiNew, Phi, FcX, FcY, FdX, FdY, Ue, Vn, U, V, DeltaX, DeltaT, Ima,Jma,nG,CFL,tMax,nTMax,D,saveNT)
    
    t = 0
    I = 0
    mom2vel!(Ue, Vn, U, V)
    save.([Phi,Ue,Vn], ["Phi", "Ue", "Vn"], I)
    log(Ima,Jma,nG,CFL,delta,DeltaT,tMax,nTMax,D,saveNT)
    while (t<tMax && I<=nTMax)
        step!(PhiNew, Phi, FcX, FcY, FdX, FdY, Ue, Vn, U, V, DeltaX, DeltaT)

        Phi = PhiNew
        t += DeltaT
        I += 1

        log(I,DeltaT,t)
        if I % saveNT == 0
            save.([Phi,Ue,Vn], ["Phi", "Ue", "Vn"], I)
        end
    end
    log()

end

simulate!(PhiNew, Phi, FcX, FcY, FdX, FdY, Ue, Vn, U, V, delta, DeltaT, Ima,Jma,nG,CFL,tMax,nTMax,D,saveNT)
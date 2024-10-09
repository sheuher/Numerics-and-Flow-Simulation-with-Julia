import Pkg; Pkg.activate("..")
using LinearAlgebra, Plots
include("mom2vel.jl")
include("calcFluxConCDS.jl")
include("calcFluxDif.jl") 
include("divMom2D.jl") 
include("poissonSolver2D.jl")
include("corrMom2D.jl") 
include("case.jl")
include("applyBCU.jl")
include("applyBCV.jl")
include("save.jl")
include("log.jl")

function step!(rhou,Ue,Vn,rhouP,divPred,p,FcXX,FcXY,FcYX,FcYY,FdXX,FdXY,FdYX,FdYY,DeltaX, DeltaT)
    nG = 1
    ImaAll, JmaAll = size(p)
	Ifi = nG+1
	Ila = ImaAll - nG
	Jfi = nG+1
	Jla = JmaAll - nG
    
    mom2vel!(Ue,Vn,rhou[:,:,1],rhou[:,:,2])

    calcFluxConCDS!(FcXX, FcXY, rhou[:,:,1], Ue, Vn, DeltaX)
    calcFluxDif!(FdXX, FdXY, rhou[:,:,1], D, DeltaX)

    rhouP[Ifi:Ila, Jfi:Jla, 1] .= rhou[Ifi:Ila, Jfi:Jla, 1] .+ DeltaT/DeltaX^2 *(
         -FcXX[Ifi:Ila, Jfi:Jla] .- FcXY[Ifi:Ila, Jfi:Jla] .+ 
          FdXX[Ifi:Ila, Jfi:Jla] .+ FdXY[Ifi:Ila, Jfi:Jla] ) 

    calcFluxConCDS!(FcYX, FcYY, rhou[:,:,2], Ue, Vn, DeltaX)
    calcFluxDif!(FdYX, FdYY, rhou[:,:,2], D, DeltaX)

    rhouP[Ifi:Ila, Jfi:Jla, 2] .= rhou[Ifi:Ila, Jfi:Jla, 2] .+ DeltaT/DeltaX^2 *(
        -FcYX[Ifi:Ila, Jfi:Jla] .- FcYY[Ifi:Ila, Jfi:Jla] .+
         FdYX[Ifi:Ila, Jfi:Jla] .+ FdYY[Ifi:Ila, Jfi:Jla] )

    divMom2D!(divPred, rhouP, DeltaX)
    poissonSolver2D!(p, divPred, DeltaT, delta, 2000, 1.0E-6)
    corrMom2D!(rhou, rhouP, p, DeltaT, delta)
    
    applyBCU!(rhou[:,:,1])
    applyBCV!(rhou[:,:,2])
end

function simulate!(rhou,Ue,Vn,rhouP,divPred,p,
                FcXX,FcXY,FcYX,FcYY,FdXX,FdXY,FdYX,FdYY,
                DeltaX, DeltaT, Ima,Jma,nG,tMax,nTMax,D,saveNT)
    
    t = 0
    I = 0
    mom2vel!(Ue, Vn, rhou[:,:,1], rhou[:,:,2])
    save.([rhou[:,:,1], rhou[:,:,2], p], ["rhou", "rhov", "p"], I)
    log(Ima,Jma,nG,NaN,delta,DeltaT,tMax,nTMax,D,saveNT)
    while (t<tMax && I<=nTMax)
        step!(rhou,Ue,Vn,rhouP,divPred,p,
                FcXX,FcXY,FcYX,FcYY,FdXX,FdXY,FdYX,FdYY,DeltaX,DeltaT)
        t += DeltaT
        I += 1

        log(I,DeltaT,t)
        if I % saveNT == 0
            save.([rhou[:,:,1], rhou[:,:,2], p], ["rhou", "rhov", "p"], I) 
        end
    end
    log()

end

simulate!(rhou,Ue,Vn,rhouP,divPred,p,
                FcXX,FcXY,FcYX,FcYY,FdXX,FdXY,FdYX,FdYY,
                delta, DeltaT, Ima,Jma,nG,tMax,nTMax,D,saveNT)
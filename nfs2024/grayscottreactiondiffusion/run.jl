import Pkg; Pkg.activate("..")
using LinearAlgebra, Plots

include("case4.jl")
include("applyBC.jl")
include("save.jl")
include("log.jl")


macro qx(ix, iy, ik) esc(:( (C[$ix+1, $iy, $ik] - C[$ix, $iy, $ik])  *D[$ik] /dx^2 )) end
macro qy(ix, iy, ik) esc(:( (C[$ix, $iy+1, $ik] - C[$ix, $iy, $ik])  *D[$ik] /dy^2 )) end

macro react(ix, iy, ik) esc(:( 
                                (-C[$ix, $iy, 1]* C[$ix, $iy, 2] ^2) * ($ik == 1) +
                                ( C[$ix, $iy, 1]* C[$ix, $iy, 2] ^2) * ($ik == 2)
                        ))end

macro kll(ix, iy, ik) esc(:( 
                                (F *(1 -C[$ix, $iy, 1])) * ($ik == 1) +
                                (-(F + k) *C[$ix, $iy, 2]) * ($ik == 2)
                        ))end

function step!(C2, C, DeltaX, DeltaT)
    nG = 1
    ImaAll, JmaAll, nPhi = size(C)
	Ifi = nG+1;     Ila = ImaAll - nG
	Jfi = nG+1;     Jla = JmaAll - nG

    dx = dy = DeltaX
    dt = DeltaT

    D = [2.0E-05, 1.0E-05]
    F = 0.029
    k = 0.057

    for k in 1:nPhi
        for j in Jfi:Jla
            for i in Ifi:Ila
                @inbounds C2[i, j, k] = C[i, j, k] + dt * ( 
                    #@react(i, j ,k) + @kll(i, j, k)
                    + @qx(i, j, k) - @qx(i-1, j, k) + @qy(i, j, k) - @qy(i, j-1, k)
                    #+ (C[i-1, j, k] + C[i+1, j, k] + C[i, j-1, k] + C[i, j+1, k] - 4*C[i, j, k]) * D[k] / DeltaX^2
                )
            end
        end
    end
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

function simulate!(PhiNew, Phi, DeltaX, DeltaT, Ima,Jma,nG,tMax,nTMax,D,saveNT)
    
    t = 0
    I = 0
    save(Phi[:,:,1], "Phi", I)
    #log(Ima,Jma,nG,CFL,delta,DeltaT,tMax,nTMax,D,saveNT)
    applyBC!(Phi)
    while (t<tMax && I<=nTMax)
        applyBC!(Phi)
        step!(PhiNew, Phi, DeltaX, DeltaT)
        Phi, PhiNew = PhiNew, Phi
        t += DeltaT
        I += 1

        log(I,DeltaT,t)
        ( I % saveNT == 0 ) &&
            #save.([Phi,Ue,Vn], ["Phi", "Ue", "Vn"], I)
            save(Phi[:,:,1], "Phi", I) 
    end
    log()

end

simulate!(PhiNew, Phi, delta, DeltaT, Ima,Jma,nG,tMax,nTMax,D,saveNT)
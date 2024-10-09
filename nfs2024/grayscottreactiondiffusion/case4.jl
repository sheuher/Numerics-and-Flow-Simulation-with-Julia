Ima = 256; Jma = 256
delta = 2.5/Ima; nG = 1

Ifi = nG +1; Ifip=Ifi+1; Ifim=Ifi-1
Ila = Ima +nG; Ilap=Ila+1; Ilam=Ila-1
ImaAll = Ima + 2*nG
Jfi = nG +1; Jfip=Ifi+1; Jfim=Ifi-1
Jla = Jma +nG; Jlap=Jla+1; Jlam=Jla-1
JmaAll = Jma + 2*nG

nPhi = 2

D = [2.0E-05, 1.0E-05]
F = 0.029
k = 0.057
tMax = 10_00_000
nTMax = 200_000
saveNT = 10_000

x = LinRange( ((-ImaAll/2) +0.5) *delta, ((ImaAll/2) -0.5) *delta, ImaAll )
y = LinRange( ((-JmaAll/2) +0.5) *delta, ((JmaAll/2) -0.5) *delta, JmaAll )
Y = ones(length(x)) * y'
X = x * ones(length(y))'

DeltaT = 0.8 * delta^2 / 2*(maximum(D)) 
#delta^2 / maximum(D) / 8.2 

# initialize scalar field
Phi  = Array{Float64}(undef, ImaAll, JmaAll, nPhi)
fill!(Phi[:,:,1], 1.0)
fill!(Phi[:,:,2], 0.0)

# 
# for j in Jfi:Jla
#     for i in Ifi:Ila
#         if ( Ima/2 - 10 < i < Ima/2 + 10 ) && ( Jma/2 - 10 < j < Jma/2 + 10 )
#             Phi[i, j, 1] = 0.5 
#             Phi[i, j, 2] = 0.25
#         end
#     end
# end
# Phi[:,:,:] += (2*rand(ImaAll, JmaAll, nPhi) .-1)/100

# RANDOM IC 
Phi[:,:,1] = ones(ImaAll, JmaAll) / 2 + (2*rand(ImaAll, JmaAll) .-1)/200
Phi[:,:,2] = ones(ImaAll, JmaAll) / 4 + (2*rand(ImaAll, JmaAll) .-1)/400

# TREFETHEN ic
Phi[:,:,1] = @. 1 - exp(-80*((x+0.05)^2 + (y'+0.05)^2))
Phi[:,:,2] = @. exp(-80*((x-0.05)^2 + (y'-0.05)^2))

# Phi[:,:,1] .= @. exp(-x^2 - y'^2) 
# Phi[:,:,2] .= 1. .- 2*Phi[:,:,1]

Phie = similar(Phi)
Phin = similar(Phi)
FdX   = similar(Phi)
FdY   = similar(Phi)
PhiNew = similar(Phi)

fill!(PhiNew, 0.0)
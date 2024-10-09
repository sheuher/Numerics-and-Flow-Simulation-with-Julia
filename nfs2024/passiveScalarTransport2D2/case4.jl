Ima = 50; Jma = 50
delta = 1/Ima; nG = 1

Ifi = nG +1; Ifip=Ifi+1; Ifim=Ifi-1
Ila = Ima +nG; Ilap=Ila+1; Ilam=Ila-1
ImaAll = Ima + 2*nG
Jfi = nG +1; Jfip=Ifi+1; Jfim=Ifi-1
Jla = Jma +nG; Jlap=Jla+1; Jlam=Jla-1
JmaAll = Jma + 2*nG

nPhi = 2

D = [1.0E-01, 1.0E-1]
F = 0.055
k = 0.062

tMax = 10.
nTMax = 1000
saveNT = 100

# velocity field: Vortex 
SF = 10.
x = LinRange( ((-ImaAll/2) +0.5) *delta, ((ImaAll/2) -0.5) *delta, ImaAll )
y = LinRange( ((-JmaAll/2) +0.5) *delta, ((JmaAll/2) -0.5) *delta, JmaAll )
Y = ones(length(x)) * y'
X = x * ones(length(y))'
U = similar(X); V = similar(V)
fill!(U, 0.0); fill!(V, 0.0)

DeltaT = delta^2 / D / 8.2
# DeltaT = CFL *delta /(maximum(U))

# initialize scalar field
Phi  = Matrix{Float64}(undef, ImaAll, JmaAll, nPhi)
nPhi = @. exp(-x^2 - y^2)
#Phi = zeros(ImaAll, JmaAll)
# for i in Ifi:Ila
#     for j in Jfi:Jla
#         Phi[i,j] = sin( (j-nG-0.5)/Ima*10*pi )
#         #Phi[i,j] = sin((i-nG-0.5)/Ima*2.0*pi) * sin((j-nG-0.5)/Jma*2.0*pi)
#     end
# end

Phie = similar(Phi)
Phin = similar(Phi)
Ue   = similar(Phi)
Vn   = similar(Phi)
FcX   = similar(Phi)
FcY   = similar(Phi)
FdX   = similar(Phi)
FdY   = similar(Phi)
PhiNew = similar(Phi)

fill!(PhiNew, 0.0)
fill!(Ue, 0.0)
fill!(Vn, 0.0)

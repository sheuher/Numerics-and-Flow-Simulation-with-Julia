Ima = 100; Jma = 100
delta = 1/Ima; nG = 1

Ifi = nG +1; Ifip=Ifi+1; Ifim=Ifi-1
Ila = Ima +nG; Ilap=Ila+1; Ilam=Ila-1
ImaAll = Ima + 2*nG
Jfi = nG +1; Jfip=Ifi+1; Jfim=Ifi-1
Jla = Jma +nG; Jlap=Jla+1; Jlam=Jla-1
JmaAll = Jma + 2*nG

DS = 1/4
CFL = 0.1
tMax = 0.1
nTMax = 10000
D = 1e-5
saveNT = 100

# velocity field: Vortex 
SF = 10.
x = LinRange( ((-ImaAll/2) +0.5) *delta, ((ImaAll/2) -0.5) *delta, ImaAll )
y = LinRange( ((-JmaAll/2) +0.5) *delta, ((JmaAll/2) -0.5) *delta, JmaAll )
Y = ones(length(x)) * y'
X = x * ones(length(y))'
R = sqrt.(X.^2 + Y.^2)
MV = SF./R
U = MV .* (Y ./R)
V = MV .* (-X ./R)

DeltaT = CFL *delta /(maximum(U))

# initialize scalar field
Phi = Matrix{Float64}(undef, ImaAll, JmaAll)
#Phi = zeros(ImaAll, JmaAll)
for i in Ifi:Ila
for j in Jfi:Jla
        Phi[i,j] = sin((i-nG-0.5)/Ima*2.0*pi) * sin((j-nG-0.5)/Jma*2.0*pi)
    end
end

Phie = similar(Phi)
Phin = similar(Phi)
Ue   = similar(Phi)
Vn   = similar(Phi)
FcX   = similar(Phi)
FcY   = similar(Phi)
FdX   = similar(Phi)
FdY   = similar(Phi)
PhiNew = similar(Phi)

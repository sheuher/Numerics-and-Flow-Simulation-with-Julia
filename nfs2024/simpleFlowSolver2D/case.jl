Ima = 15; Jma = 35
delta = 1/Ima; nG = 1

Ifi = nG +1; Ifip=Ifi+1; Ifim=Ifi-1
Ila = Ima +nG; Ilap=Ila+1; Ilam=Ila-1
ImaAll = Ima + 2*nG
Jfi = nG +1; Jfip=Ifi+1; Jfim=Ifi-1
Jla = Jma +nG; Jlap=Jla+1; Jlam=Jla-1
JmaAll = Jma + 2*nG

#CFL = 0.1
DeltaT = 0.001
tMax = 1.0
nTMax = 500
D = 1.0E-1
saveNT = 5

# initialize momentom field
rhou = zeros(ImaAll, JmaAll, 2)
rhou[:,:,1] .= 1

# initialize pressure field
p = zeros(ImaAll, JmaAll)

# other fields initiliazation
rhouP = zeros(size(rhou)...)
divPred = zeros(size(p)...)
Ue   = zeros(size(p)...)
Vn   = zeros(size(p)...)
FcXX   = zeros(size(p)...)
FcYX   = zeros(size(p)...)
FdXX   = zeros(size(p)...)
FdYX   = zeros(size(p)...)
FcXY   = zeros(size(p)...)
FcYY   = zeros(size(p)...)
FdXY   = zeros(size(p)...)
FdYY   = zeros(size(p)...)
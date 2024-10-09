function [X,resisdual,nitr] = PoissonCG(f, BCtype)
% solve pressure poisson equation using Conjugate Gradient method Laplacian(P) = f
% PE + PW + PN + PS -4PC = rhs*Delta^2
% B.C. zero gradient at all boundaries and 1 point known to be zero i.e. P(Ifi,Jfi) = 0


global ImaAll Ifi Ila Ifim Ilap 
global JmaAll Jfi Jla Jfim Jlap 

%solver parameters
tol = 1e-6;
maxitr= 2000;

%diag factors deals with B.Cs. 
diagFactor = diagFactor_PP(BCtype); 
%invdiagFactor = -1*diagFactor.^-1; % (-1)--> diag factors are positive while in the equation are negative. 

X = zeros(ImaAll,JmaAll); %x = zeros(length(b),1);


                                    
r = res(X,f,diagFactor); % r(ImaAll,JmaAll)            %r = b - A * x;  
p = r;                 % p(ImaAll,JmaAll)              %p = r;  
rsold = sum(r.*r,'all');                               %rsold = r' * r; 

for i = 1:maxitr
    Ap = leftMultA(p,diagFactor);                      %Ap = A * p; 
    alpha = rsold / (sum(Ap.*p,'all')); 
    X = X + alpha * p;                                 
    r = r - alpha * Ap;                                
    rsnew = sum(r.*r,'all');                           
    if sqrt(rsnew) < tol
        break
    end
    p = r + (rsnew / rsold) * p;                       
    rsold = rsnew; 
    
    nitr = i;
end

%set zero gradient B.C
X(Ifim,:) = X(Ifi,:);
X(Ilap,:) = X(Ila,:);
X(:,Jfim) = X(:,Jfi);
X(:,Jlap) = X(:,Jla);

resisdual = sqrt(rsnew);


end



function RES = res(X,rhs,diagFactor)
%calculate residulas of the solver
global ImaAll Ifi Ila Ifim Ilap Ifip Ilam
global JmaAll Jfi Jla Jfim Jlap Jfip Jlam
global Delta
RES = zeros(ImaAll,JmaAll);
Ax = leftMultA(X,diagFactor);
RES(Ifi:Ila,Jfi:Jla) = Delta^2* rhs(Ifi:Ila,Jfi:Jla) - Ax(Ifi:Ila,Jfi:Jla);

end

function Ax = leftMultA(X,diagFactor)
%left multipication of matrix A and vector P: A*p in matrix form
% A(n*n,n*n) , P(n,n)

global ImaAll Ifi Ila Ifim Ilap Ifip Ilam
global JmaAll Jfi Jla Jfim Jlap Jfip Jlam
global Delta
Ax =zeros(ImaAll,JmaAll);

Ax(Ifi:Ila,Jfi:Jla) = X(Ifim:Ilam,Jfi:Jla) + X(Ifip:Ilap,Jfi:Jla)...
    + X(Ifi:Ila,Jfim:Jlam) + X(Ifi:Ila,Jfip:Jlap)...
    -diagFactor(Ifi:Ila,Jfi:Jla).*X(Ifi:Ila,Jfi:Jla);

end


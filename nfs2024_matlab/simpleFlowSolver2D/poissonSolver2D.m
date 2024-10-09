% # net.m ###############################################
% Programm: 		    	simpleFlowSolver2D
% Content:					poisson solver
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% function that solve the poisson equation for pressure
% correction
% =======================================================

function [p, nIt, eps] = poissonSolver2D(p, divPred, ...
    A, BCs) % ###########################################
% ------- parameters ------------------------------------
nItMax = 10000;
epsMax = 1.0e-6;
nIt = 0;
eps = 10.0;
global Ifi Ifim Ifip Ila Ilam Ilap 
global Jfi Jfim Jfip Jla Jlam Jlap
global Delta DeltaT ImaAll JmaAll
% ------- preallocating and preindexing -----------------
RHS = Delta^2 /DeltaT *divPred(Ifi:Ila, Jfi:Jla);
diagFac = ones(ImaAll,JmaAll)*4;
diagFac(Ifi,:) = diagFac(Ifi,:)-1; 
diagFac(:,Jfi) = diagFac(:,Jfi)-1;
diagFac(:,Jla) = diagFac(:,Jla)-1;
diagFac= 1./diagFac;
% ------- operation jacobi method -----------------------
while( nIt<=nItMax && eps>epsMax )
    pold = p;
    % please note that:
    % pE = p(Ifip:Ilap, Jfi:Jla)
    % pW = p(Ifim:Ilam, Jfi:Jla)
    % pN = p(Ifi:Ila, Jfip:Jlap)
    % pS = p(Ifi:Ila, Jfim:Jlam)
    % pC = p(Ifi:Ila, Jfi:Jla)
    p(Ifi:Ila, Jfi:Jla) = diagFac(Ifi:Ila, Jfi:Jla) .* (...
        pold(Ifip:Ilap, Jfi:Jla) + ...
        pold(Ifim:Ilam, Jfi:Jla) + ...
        pold(Ifi:Ila, Jfip:Jlap) +...
        pold(Ifi:Ila, Jfim:Jlam) - RHS);
    % BC! for pressure
    %#example zero gradient at the inlet
    %# fixed pressure of 0 at 3 other outlet
    % homogeneus Neumann BC everywhere except top
    if BCs(1) == 0
        p(Ilap, Jfi:Jla) = 0;
    else 
        p(Ilap, Jfi:Jla) = p(Ila, Jfi:Jla);
    end
    if BCs(2) == 0
        p(Ifi:Ila, Jfim) = 0;
    else
        p(Ifi:Ila, Jfim) = p(Ifi:Ila, Jfi);
    end
    if BCs(3) == 0
        p(Ifim, Jfi:Jla) = 0;
    else
        p(Ifim, Jfi:Jla) = p(Ifi, Jfi:Jla);
    end
    if BCs(4) == 0
        p(Ifi:Ila, Jlap) = 0;
    else
        p(Ifi:Jla, Jlap) = p(Ifi:Jla, Jla);
    end
    % update jacobi step
    nIt = nIt +1;
    eps = max( abs(p - pold),[],'all' );
end

 


end % ###################################################

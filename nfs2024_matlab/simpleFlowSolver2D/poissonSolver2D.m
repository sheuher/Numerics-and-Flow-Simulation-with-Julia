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

function [p, nIt, eps] = poissonSolver2D(p, divPred) % ##
% ------- parameters ------------------------------------
nItMax = 10000;
epsMax = 1.0e-5;
nIt = 0;
eps = 10.0;
global Ifi Ifim Ifip Ila Ilam Ilap 
global Jfi Jfim Jfip Jla Jlam Jlap
global Delta DeltaT
% ------- preallocating and preindexing -----------------
RHS = Delta^2 /DeltaT *divPred(Ifi:Ila, Jfi:Jla);
% ------- operation jacobi method -----------------------
while( nIt<=nItMax && eps>epsMax )
    pold = p;
    % please note that:
    % pE = p(Ifip:Ilap, Jfi:Jla)
    % pW = p(Ifim:Ilam, Jfi:Jla)
    % pN = p(Ifi:Ila, Jfip:Jlap)
    % pS = p(Ifi:Ila, Jfim:Jlam)
    % pC = p(Ifi:Ila, Jfi:Jla)
    p(Ifi:Ila, Jfi:Jla) = 0.25 * (...
        pold(Ifip:Ilap, Jfi:Jla) + ...
        pold(Ifim:Ilam, Jfi:Jla) + ...
        pold(Ifi:Ila, Jfip:Jlap) +...
        pold(Ifi:Ila, Jfim:Jlam) - RHS);
    % BC! for pressure
    %#example zero gradient at the inlet
    %# fixed pressure of 0 at 3 other outlet
    % homogeneus Neumann BC everywhere except top
    p(Ifim, Jfi:Jla) = 0;%p(Ifi, Ifi:Ila);
    p(Ilap, Jfi:Jla) = 0;%p(Ila, Jfi:Jla);
    p(Ifi:Ila, Jfim) = 0;%p(Ifi:Ila, Jfi);
    p(Ifi:Ila, Jlap) = 0;
    % update jacobi step
    nIt = nIt +1;
    eps = max( abs(p - pold),[],'all' );
end
end % ###################################################
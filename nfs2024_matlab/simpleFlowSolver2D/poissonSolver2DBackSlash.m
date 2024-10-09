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

function [p, nIt, eps] = poissonSolver2DBackSlash(p, ...
    divPred, A, BCs) % ##################################
% ------- parameters ------------------------------------
global Ifi Ila 
global Jfi Jla 
global Delta DeltaT
global Jma Ima
% ------- preallocating and preindexing -----------------
RHS = Delta^2 /DeltaT *divPred(Ifi:Ila, Jfi:Jla);
% ----- solve with x=A\b --------------------------------
M = diag(1./diag(A));
RHS = reshape(RHS, [Ima*Jma,1]);
% [pp, ~, eps, nIt] = pcg(A,RHS,1e-6,10000);
p(Ifi:Ila,Jfi:Jla) = reshape(A\RHS, [Ima,Jma]);
eps=0; nIt=-1;
end % ###################################################

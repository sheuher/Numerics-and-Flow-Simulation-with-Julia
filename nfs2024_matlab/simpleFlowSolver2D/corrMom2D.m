% # net.m ###############################################
% Programm: 		    	simpleFlowSolver2D
% Content:					apply p corr to Mom field
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% function to apply momentum field correction with 
% pressure gradient
% =======================================================

function rhou = corrMom2D(rhou, rhouP, p)% ############
% ------- parameters ------------------------------------
global Ifi Ifim Ifip Ila Ilam Ilap 
global Jfi Jfim Jfip Jla Jlam Jlap
global Delta DeltaT
% ------- preallocating and preindexing -----------------
pE = p(Ifip:Ilap, Jfi:Jla);
pW = p(Ifim:Ilam, Jfi:Jla);
pN = p(Ifi:Ila, Jfip:Jlap);
pS = p(Ifi:Ila, Jfim:Jlam);
rhouPX = rhouP(Ifi:Ila, Jfi:Jla, 1);
rhouPY = rhouP(Ifi:Ila, Jfi:Jla, 2);
% ------- operation -------------------------------------
rhou(Ifi:Ila, Jfi:Jla, 1) = rhouPX - DeltaT *0.5 ...
    /Delta * ( pE - pW );
rhou(Ifi:Ila, Jfi:Jla, 2) = rhouPY - DeltaT *0.5 ...
    /Delta * ( pN - pS );
end % ###################################################

% # net.m ###############################################
% Programm: 		    	simpleFlowSolver2D
% Content:					calculate div of momentum
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% function that calculate the divergence of the momentum
% fields before corrected
% =======================================================

function divPred = divMom2D(divPred, rhouP) % ###########
% ------- parameters ------------------------------------
global Ifi Ifim Ifip Ila Ilam Ilap 
global Jfi Jfim Jfip Jla Jlam Jlap
global Delta
% ------- preallocating and preindexing -----------------
rhouPE = rhouP(Ifip:Ilap, Jfi:Jla, 1);
rhouPW = rhouP(Ifim:Ilam, Jfi:Jla, 1);
rhouPN = rhouP(Ifi:Ila, Jfip:Jlap, 2);
rhouPS = rhouP(Ifi:Ila, Jfim:Jlam, 2);
% ------- operation -------------------------------------
divPred(Ifi:Ila, Jfi:Jla) = 0.5 /Delta * ( ...
    rhouPE - rhouPW + rhouPN - rhouPS );
end % ###################################################
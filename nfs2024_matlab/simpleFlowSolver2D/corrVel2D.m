% # net.m ###############################################
% Programm: 		    	simpleFlowSolver2D
% Content:					apply p corr to Velo field
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% function to apply velocity field correction from 
% pressure gradient on cell surfaces
% =======================================================

function [Ue, Vn] = corrVel2D(Ue, Vn, p)% ###############
% ------- parameters ------------------------------------
rho = 1                                                 ;
global Ifi Ifim Ifip Ila Ilam Ilap 
global Jfi Jfim Jfip Jla Jlam Jlap
global Delta DeltaT
% ------- preallocating and preindexing -----------------
pCX = p(Ifi:Ilap, Jfi:Jla);
pW =  p(Ifim:Ila, Jfi:Jla);
pCY = p(Ifi:Ila, Jfi:Jlap);
pS =  p(Ifi:Ila, Jfim:Jla);
% ------- operation -------------------------------------
Ue(Ifim:Ila, Jfi:Jla) = Ue(Ifim:Ila, Jfi:Jla) - ...
    DeltaT/rho/Delta * ( pCX - pW );
Vn(Ifi:Ila, Jfim:Jla) = Vn(Ifi:Ila, Jfim:Jla) - ...
    DeltaT/rho/Delta * ( pCY - pS );
end % ###################################################
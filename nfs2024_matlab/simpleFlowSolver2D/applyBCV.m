% # applyBCV.m ##########################################
% Programm: 		    	simpleFlowSolver2D
% Content:					apply boundary conditions
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% functions for boundary condition applications
%    - application of periodic boundary condition to
%      to both x, y directions
% =======================================================

function Phi = applyBCV(Phi) % ##########################
% ------- parameters ------------------------------------
global Ifi Ifim Ila Ilap Jfi Jfim Jla Jlap
% ------- operattions -----------------------------------
% Dirichlet on all sides
Phi(Ifim, Jfi:Jla) = 0;
Phi(Ilap, Jfi:Jla) = 0;
Phi(Ifi:Ila, Jfim) = 0;
Phi(Ifi:Ila, Jlap) = 0;
end % ###################################################
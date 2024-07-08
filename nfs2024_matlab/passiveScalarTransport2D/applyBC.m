% # applyBC.m ###########################################
% Programm: 		    	passiveScalarTranport2D
% Content:					apply boundary conditions
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% functions for boundary condition applications
%    - application of periodic boundary condition to
%      to both x, y directions
% =======================================================

function Phi = applyBC(Phi) % ###########################
% ------- parameters ------------------------------------
global Ifi Ifim Ila Ilap Jfi Jfim Jla Jlap
% ------- operattions -----------------------------------
Phi(Ifim, Jfi:Jla) = Phi(Ila, Jfi:Jla);
Phi(Ilap, Jfi:Jla) = Phi(Ifi, Jfi:Jla);
Phi(Ifi:Ila, Jfim) = Phi(Ifi:Ila, Jla);
Phi(Ifi:Ila, Jlap) = Phi(Ifi:Ila, Jfi);
end % ###################################################
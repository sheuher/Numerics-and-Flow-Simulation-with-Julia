% # net.m ###############################################
% Programm: 		    	simpleFlowSolver2D
% Content:					domain size and indexing 
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% initiate indexing system for specified domain size
% =======================================================

% #######################################################
% ------- parameters: initiate globally -----------------
global Ima Jma nG
global Ifi Ifim Ifip
global Ila Ilam Ilap Jfi Jfim Jfip
global Jla Jlam Jlap
global ImaAll JmaAll
% ------- parameters: define grid -----------------------
Ima =    25                                             ;
Jma =    25                                            ;
nG =     1                                              ;
% ------- parameters: indexing for grid -----------------
Ifi = nG +1;   Ifim = Ifi -1; Ifip = Ifi +1;
Ila = Ima +nG; Ilam = Ila -1; Ilap = Ila +1;
Jfi = nG +1;   Jfim = Jfi -1; Jfip = Jfi +1;
Jla = Jma +nG; Jlam = Jla -1; Jlap = Jla +1;
ImaAll = Ima + 2*nG;
JmaAll = Jma + 2*nG;
% #######################################################
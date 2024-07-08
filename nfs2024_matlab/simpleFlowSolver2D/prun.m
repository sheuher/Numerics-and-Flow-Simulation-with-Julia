% # run.m ###############################################
% Programm: 		    	simpleFlowSolver2D
% Content:					define simulation params
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% define simulation params such as CFL, Dx, Dt, tmax
% =======================================================

% #######################################################
global Delta DeltaT D saveNT nTMax tMax SAVE
% ------- parameters: run setup -------------------------
Delta =     1.0/Ima                                     ;
CFL =       0.1                                         ;
DeltaT =    0.001 ;% xmax(abs([U;V]),[],'all');
D =         1.0E-1                                      ;
saveNT =    10                                          ;
nTMax =     500                                        ;
tMax =      100.0                                       ;
SAVE =      false                                       ;
% #######################################################

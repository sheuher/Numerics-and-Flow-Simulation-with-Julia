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
DeltaT =    CFL*Delta ;% xmax(abs([U;V]),[],'all');
D =         1.0E-2                                      ;
saveNT =    10                                          ;
nTMax =     1000                                         ;
tMax =      100.0                                       ;
SAVE =      false                                       ;
% #######################################################
% 1/dt < D/delta^2 --> dt < delta**2/D
% dt < dx/Umax

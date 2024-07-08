% # run.m ###############################################
% Programm: 		    	passiveScalarTranport2D
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
Delta =     1/Ima                                       ;
CFL =       0.1                                         ;
DeltaT = CFL *Delta /max([U;V],[],'all');
D =         1.0E-5                                      ;
saveNT =    10                                          ;
nTMax =     10000                                       ;
tMax =      10.0                                        ;
SAVE =      false                                       ;
% #######################################################
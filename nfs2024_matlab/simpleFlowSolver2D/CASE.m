% # CASE.m ##############################################
% Programm: 		    	simpleFlowSolver2D
% Content:					define case, initiate fields
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% initiate fields with the defined case
% =======================================================

% #######################################################
rho  =  1                                               ;
dims =  [ImaAll, JmaAll];
% ------- parameters: mom & p fields initialization -----
rhou = zeros([dims 2]);
rhou(:, :, 1) = 1;
p =     zeros(dims);
BCs = [1,1,1,0];
% ------- parameters: other fields initialization -------
rhouP = zeros([dims 2]);
divPred=zeros(dims);
Ue =    zeros(dims);
Vn =    zeros(dims);
FcXX =  zeros(dims);
FcYX =  zeros(dims);
FdXX =  zeros(dims);
FdYX =  zeros(dims);
FcXY =  zeros(dims);
FcYY =  zeros(dims);
FdXY =  zeros(dims);
FdYY =  zeros(dims);
A    =  zeros(Ima*Jma, Ima*Jma);
% #######################################################
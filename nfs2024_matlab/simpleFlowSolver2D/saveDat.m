% # net.m ###############################################
% Programm: 		    	simpleFlowSolver2D
% Content:					utility tool save field
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% save field data to .dat file for analysis
% =======================================================

function r = saveDat(Phi, fileName) % ###################
% ------- parameters ------------------------------------
global SAVE
filePath = sprintf("dat/%s.dat", fileName);
% ------- operation -------------------------------------
if(SAVE == false)
    return
end
save(filePath, 'Phi', '-ascii');
r = 1;
end
% #######################################################

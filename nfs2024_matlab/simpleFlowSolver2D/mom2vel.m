% # net.m ###############################################
% Programm: 		    	simpleFlowSolver2D
% Content:					velocity field interpolation
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% function to interpolate velocity field u[:,:], and 
% v[:,:] on CELL SURFACES from the momentum fields 
% rhou[:,:] and rhov[:,:] stored on cell centre
% =======================================================

function [u, v] = mom2vel(u,v, rhou, rhov) % ############
% ------- parameters ------------------------------------
rho = 1;
global Ifi Ifim Ila Ilap Jfi Jfim Jla Jlap
% ------- preallocating and preindexing -----------------
rhouE = rhou(Ifi:Ilap, :);
rhouW = rhou(Ifim:Ila, :);
rhovN = rhov(:, Jfi:Jlap);
rhovS = rhov(:, Jfim:Jla);
% ------- operation -------------------------------------
u(Ifim:Ila, :) = 1/rho *0.5 * (rhouE + rhouW);
v(:, Jfim:Jla) = 1/rho *0.5 * (rhovN + rhovS);
end % ###################################################
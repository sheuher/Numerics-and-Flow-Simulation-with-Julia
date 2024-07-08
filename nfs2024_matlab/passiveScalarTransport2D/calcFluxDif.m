% # net.m ###############################################
% Programm: 		    	passiveScalarTranport2D
% Content:					diffusive fluxes clac CDS 
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% function for diffusive fluxes calculation with CDS
% differencing scheme
% =======================================================

function [fluxDifX, fluxDifY] = calcFluxDif( ...
    fluxDifX,fluxDifY, Phi) % ###########################
% ------- parameters ------------------------------------
global Ifi Ifim Ifip Ila Ilam Ilap
global Jfi Jfim Jfip Jla Jlam Jlap
global Delta D
% ------- preallocating and preindexing -----------------
PhiC = Phi(Ifi:Ila, Jfi:Jla);
PhiE = Phi(Ifip:Ilap, Jfi:Jla);
PhiW = Phi(Ifim:Ilam, Jfi:Jla);
PhiN = Phi(Ifi:Ila, Jfip:Jlap);
PhiS = Phi(Ifi:Ila, Jfim:Jlam);
% ------- operation -------------------------------------
fluxDifX(Ifi:Ila, Jfi:Jla) = D *Delta * ( ...
    PhiE + PhiW - 2*PhiC);
fluxDifY(Ifi:Ila, Jfi:Jla) = D *Delta * ( ...
    PhiN + PhiS - 2*PhiC);
% #######################################################

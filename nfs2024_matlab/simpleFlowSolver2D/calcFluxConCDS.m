% # net.m ###############################################
% Programm: 		    	simpleFlowSolver2D
% Content:					convective fluxes clac CDS 
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% function for convective fluxes calculation with CDS
% differencing scheme:
%   - Phi(:,:) is stored in face centre; 
%   - U(:,:) and V(:,:) are stored in face surfaces
% =======================================================

function [fluxConX, fluxConY] = calcFluxConCDS( ... 
    fluxConX,fluxConY, Phi, U, V) % #####################
% ------- parameters ------------------------------------
global Ifi Ifim Ifip Ila Ilam Ilap
global Jfi Jfim Jfip Jla Jlam Jlap
global Delta
% ------- preallocating and preindexing -----------------
PhiC = Phi(Ifi:Ila, Jfi:Jla);
PhiE = Phi(Ifip:Ilap, Jfi:Jla);
PhiW = Phi(Ifim:Ilam, Jfi:Jla);
PhiN = Phi(Ifi:Ila, Jfip:Jlap);
PhiS = Phi(Ifi:Ila, Jfim:Jlam);

Ue = U(Ifi:Ila, Jfi:Jla);
Uw = U(Ifim:Ilam, Jfi:Jla);
Vn = V(Ifi:Ila, Jfi:Jla);
Vs = V(Ifi:Ila, Jfim:Jlam);
% ------- operation -------------------------------------
fluxConX(Ifi:Ila, Jfi:Jla) = 0.5 *Delta * ( ...
    (PhiE +PhiC).*Ue - (PhiC +PhiW).*Uw );
fluxConY(Ifi:Ila, Jfi:Jla) = 0.5 *Delta * ( ...
    (PhiN +PhiC).*Vn - (PhiC +PhiS).*Vs );
end % ###################################################

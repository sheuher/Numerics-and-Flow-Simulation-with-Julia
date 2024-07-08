% # net.m ###############################################
% Programm: 		    	passiveScalarTranport2D
% Content:					iterate one timestep
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% function that perform iteration for one time step
% =======================================================

function [PhiNew] = step(PhiNew, Phi, FcX, FcY, ...
    FdX, FdY, Ue, Vn)% ##################################
% ------- parameters ------------------------------------
global Ifi Ila Jfi Jla Delta DeltaT 
% ------- operations ------------------------------------
% step
%  \_ applyBC x many
Ue = applyBC(Ue); Vn = applyBC(Vn); Phi = applyBC(Phi);
%  \_ calcFluxConCDS
[FcX, FcY] = calcFluxConCDS(FcX,FcY,Phi,Ue,Vn);
%  \_ calcFluxDif
[FdX, FdY] = calcFluxDif(FdX,FdY,Phi);
%  \_ PhiNew = Phi + ...
PhiNew(Ifi:Ila, Jfi:Jla) = Phi(Ifi:Ila, Jfi:Jla) + ...
    DeltaT/Delta^2 .* ( - FcX(Ifi:Ila, Jfi:Jla) - ...
    FcY(Ifi:Ila, Jfi:Jla) + FdX(Ifi:Ila, Jfi:Jla) + ...
    FdY(Ifi:Ila, Jfi:Jla)  );
%  \_ applyBC
PhiNew = applyBC(PhiNew);
end % ###################################################

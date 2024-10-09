% # net.m ###############################################
% Programm: 		    	simpleFlowSolver2D
% Content:					iterate one timestep
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% function that perform iteration for one time step
% =======================================================

function [rhou,Ue,Vn,p,nIt,eps] = step(rhou,Ue,Vn, ...
    rhouP, divPred, p, FcXX, FcXY, FcYX, FcYY, ...
    FdXX, FdXY, FdYX, FdYY, A, BCs) % ###################
% ------- parameters ------------------------------------
global Ifi Ila Jfi Jla Delta DeltaT 
% ------- operations ------------------------------------
% step 
%  \_ calcFluxConCDS
%  \_ calcFluxDif
%  \_ apply fluxes
%  \_ applyBC (rhouP)
%  \_ divMom2D
%  \_ poissonSolver2D
%  \_ corrMom2D
%  \_ applyBC (rhou)
%  \_ mom2vel (rhouP) & corrVel2D (update Ue,Vn)
[Ue,Vn] = mom2vel(Ue,Vn,rhou(:,:,1),rhou(:,:,2));

[FcXX,FcXY]= calcFluxConUDS(FcXX,FcXY,rhou(:,:,1),Ue,Vn);
[FdXX,FdXY]= calcFluxDif(FdXX,FdXY,rhou(:,:,1));
rhouP(Ifi:Ila,Jfi:Jla, 1) = rhou(Ifi:Ila,Jfi:Jla,1) + ...
    DeltaT/Delta^2 * ( - FcXX(Ifi:Ila,Jfi:Jla) - ...
    FcXY(Ifi:Ila,Jfi:Jla) + FdXX(Ifi:Ila,Jfi:Jla) + ...
    FdXY(Ifi:Ila,Jfi:Jla) );

[FcYX,FcYY]= calcFluxConCDS(FcYX,FcYY,rhou(:,:,2),Ue,Vn);
[FdYX,FdYY]= calcFluxDif(FdYX,FdYY,rhou(:,:,2));
rhouP(Ifi:Ila,Jfi:Jla, 2) = rhou(Ifi:Ila,Jfi:Jla,2) + ...
    DeltaT/Delta^2 * ( - FcYX(Ifi:Ila,Jfi:Jla) - ...
    FcYY(Ifi:Ila,Jfi:Jla) + FdYX(Ifi:Ila,Jfi:Jla) + ...
    FdYY(Ifi:Ila,Jfi:Jla) );

[rhouP(:,:,1)] = applyBCU(rhouP(:,:,1));
[rhouP(:,:,2)] = applyBCV(rhouP(:,:,2));

[divPred] = divMom2D(divPred,rhouP);
%[p,nIt,eps] = poissonSolver2D(p,divPred,A,BCs);
[p,nIt,eps] = poissonSolver2DBackSlash(p,divPred,A,BCs);
[rhou] = corrMom2D(rhou,rhouP,p);

[rhou(:,:,1)] = applyBCU(rhou(:,:,1));
[rhou(:,:,2)] = applyBCV(rhou(:,:,2));
end % ###################################################
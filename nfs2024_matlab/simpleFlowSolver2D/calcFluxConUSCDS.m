% # calcFluxConUSCDS.m ##################################
% Programm: 		    	passiveScalarTranport2D
% Content:					convective fluxes clac CDS 
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% PLEASE VERIFY FOR CORRECTNESS!!!!!!
% function for convective fluxes calculation with USCDS
% differencing scheme:
%   - Phi(:,:) is stored in face centre; 
%   - U(:,:) and V(:,:) are stored in face surfaces
% =======================================================

function [fluxConX, fluxConY] = calcFluxConUSCDS( ... 
    fluxConX,fluxConY, Phi, U, V) % #####################
% ------- parameters ------------------------------------
global Ifi Ifim Ifip Ila Ilam Ilap
global Jfi Jfim Jfip Jla Jlam Jlap
global Delta DeltaT Ima Jma
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

Fe = zeros(Ima, Jma);           % fluxes arrays
Fw = zeros(Ima, Jma);
Fn = zeros(Ima, Jma);
Fs = zeros(Ima, Jma);
% ------- operation -------------------------------------
Fe(:,:) = (1 + Ue *DeltaT/Delta) .*PhiC +...
    (1 - Ue *DeltaT/Delta) .*PhiE;
Fw(:,:) = (1 + Uw *DeltaT/Delta) .*PhiW +...
    (1 - Uw *DeltaT/Delta) .*PhiC;

Fn(:,:) = (1 + Vn *DeltaT/Delta) .*PhiC +...
    (1 - Vn *DeltaT/Delta) .*PhiN;
Fs(:,:) = (1 + Vs *DeltaT/Delta) .*PhiS +...
    (1 - Vs *DeltaT/Delta) .* PhiC;

fluxConX(Ifi:Ila, Jfi:Jla) = 0.5*Delta *(Fe.*Ue - Fw.*Uw);
fluxConY(Ifi:Ila, Jfi:Jla) = 0.5*Delta *(Fn.*Vn - Fs.*Vs);
end % ###################################################
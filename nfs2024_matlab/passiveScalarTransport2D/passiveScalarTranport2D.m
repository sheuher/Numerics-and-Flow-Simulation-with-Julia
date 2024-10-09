% # passiveScalarTransport2D.m ##########################
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% a beautiful scalar field transport predictor
% =======================================================

% #######################################################
% ------- define params grid and indexes ----------------
pnet     % net gives us the familiar Ifi, Ila... indexes
% ------- define case -----------------------------------
CASE    % CASE initiates all required field matrixes
% ------- define params for execution of the calculation
prun
% ------- test graphically defined velocity field -------
%quiver(X,Y,U,V,0.5)
%heatmap(V)
% ------- simulate! -------------------------------------
simulate(PhiNew, Phi, FcX, FcY, FdX, FdY, Ue, Vn, U, V);
% #######################################################

% #######################################################
% ------- functions defined -----------------------------
% [Ue, Vn] = mom2vel(Ue,Vn,U,V);
% [FcX, FcY] = calcFluxConCDS(FcX,FcY,Phi,U,V);
% [FdX, FdY] = calcFluxDif(FdX,FdY,Phi);
% [Phi] = applyBC(Phi);
% [PhiNew] = step(PhiNew,Phi,FcX,FcY,FdX,FdY,Ue,Vn);
    % step
    %  \_ applyBC x many
    %  \_ calcFluxConCDS
    %  \_ calcFluxDif
    %  \_ PhiNew = Phi + ...
    %  \_ applyBC
% #######################################################

function r = simulate(PhiNew, Phi, FcX, FcY, FdX, ...
    FdY, Ue, Vn, U, V) % ################################
        % simulate 
        %  \_ mom2vel
        %  \_ applyBC x many
        %  \_ save
        %  \_ log
        %  \_ whileloop
        %       \_ step
        %       \_ Phi = PhiNew
        %       \_ t += DeltaT
        %       \_ I += 1
        %       \_ log
        %       \_ if I % saveNT == 0: save
        %  \_ log
% ------- parameters ------------------------------------
global DeltaT saveNT nTMax tMax
t = 0.0;
I = 0;
% ------- SIMULATE! -------------------------------------
% ------- mom2vel 
[Ue, Vn] = mom2vel(Ue,Vn,U,V);
% ------- applyBCs 
Phi = applyBC(Phi); Ue = applyBC(Ue); Vn = applyBC(Vn);
% ------- save 
saveDat(Phi, 'Phi0');
% ------- log 
fprintf("SIMULATION STARTED")
% ------- whileloop 
while( t<tMax && I<=nTMax ) 
    % ------- step 
    PhiNew = step(PhiNew,Phi,FcX,FcY,FdX,FdY,Ue,Vn);
    % ------- Phi = PhiNew 
    Phi = PhiNew;
    % ------- t += DeltaT 
    t = t +DeltaT;
    % ------- I += 1 
    I = I +1;
    % ------- log 
    fprintf("I = %d, t = %.04f s\n", I, t)
    % ------- if I % saveNT == 0: save 
    if(mod(I, saveNT) == 0)
        saveDat(Phi, sprintf('Phi%i', I))
        %surf(Phi); zlim([-1,1]); colormap gray; 
        %shading faceted;
        heatmap(Phi, 'Colormap', gray, 'ColorLimits', ...
            [-1.2 1.2])
        %
         grid off
        % axis normal
        drawnow
        
    end
% ------- log 
end
r = 1;
end % ###################################################
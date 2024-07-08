% # simpleFlowSolver2D.m ################################
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% a beautiful and simple flow solver
% =======================================================

% #######################################################
% ------- define params grid and indexes ----------------
pnet      % net gives us the familiar Ifi, Ila... indexes
% ------- initialization --------------------------------
CASE      % CASE initiates all required field matrixes
% ------- define params for execution of the calculation
prun
% ------- test graphically momentum fields --------------
% quiver(X,Y,.,.,0.5)% comment out simulate below to test
% heatmap(...)
%surface(rhou(:,:,1))
% ------- SIMULATE! -------------------------------------
simulate(rhou,Ue,Vn, rhouP, divPred, p, ...
    FcXX, FcXY, FcYX, FcYY, FdXX, FdXY, FdYX, FdYY);
% #######################################################

% #######################################################
% ------- functions defined -----------------------------
% [rhou(:,:,1)] = applyBCU(rhou(:,:,1));
% [rhou(:,:,2)] = applyBCV(rhou(:,:,2));
% [Ue,Vn] =     mom2vel(Ue,Vn,rhou(:,:,1),rhou(:,:,2));
% [FcXX,FcYX] = calcFluxConCDS(FcXX,FcYX,rhou(:,:,1), ...
%                                Ue,Vn);
% [FdXX,FdYX] = calcFluxDif(FdXX,FdYX,rhou(:,:,1));
% [Phi] =       applyBC(Phi);
% [divPred] =   divMom2D(divPred,rhouP);
% [p,nIt,eps] = poissonSolver2D(p,divPred)
% [rhou] =      corrMom2D(rhou,rhouP,p)
% [Ue,Vn] =     corrVel2D(Ue,Vn,p)
% [rhou,Ue,Vn,nIt,eps] = step(rhou,Ue,Vn,rhouP, ...
%   divPred,p,FcXX,FcXY,FcYX,FcYY,FdXX, FdXY, FdYX, ...
%   FdYY);
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
% #######################################################

function r = simulate(rhou,Ue,Vn, rhouP, divPred, p, ...
    FcXX, FcXY, FcYX, FcYY, FdXX, FdXY, FdYX, FdYY)  % ##
        % simulate 
        %  \_ applyBC x many
        %  \_ mom2vel calc U,V from RhoU,RhoV
        %  \_ save
        %  \_ log start
        %  \_ whileloop
        %       \_ step
        %       \_ rhou = rhouNew       !update momentum
        %       \_ Ue,Vn = UeNew,VnNew  !update face Vel
        %       \_ t += DeltaT
        %       \_ I += 1
        %       \_ log
        %       \_ if I % saveNT == 0: save, visualise
        %  \_ log end
% ------- parameters ------------------------------------
global DeltaT saveNT nTMax tMax
global Ifi Ila Jfi Jla 
t = 0.0;
I = 0;
% ------- SIMULATE! -------------------------------------
% ------- applyBCs
[rhou(:,:,1)] = applyBCU(rhou(:,:,1));
[rhou(:,:,2)] = applyBCV(rhou(:,:,2));
% ------- mom2vel 
[Ue, Vn] = mom2vel(Ue, Vn, rhou(:,:,1), rhou(:,:,2));
% ------- save initial field
fprintf("SIMULATION STARTED\n")
while( t<tMax && I<=nTMax ) 
    % --- step, update momentum & face velocity
    [rhou,Ue,Vn,p,nIt,eps] = step(rhou,Ue,Vn,rhouP, ...
        divPred,p,FcXX,FcXY,FcYX,FcYY,FdXX, FdXY, ...
        FdYX, FdYY);
    t = t +DeltaT;
    I = I +1;
    fprintf("I = %d, t = %.04f s, ", I, t)
    fprintf("nJac = %d, eJac = %.03d\n", nIt, eps)
    if(mod(I, saveNT) == 0)
        %surface(p');
        surface(rhou(:,:,1)');
         % surface(p(:,:)');
        %heatmap(p');
        %heatmap(rhou(Ifi:Ila,Jfi:Jla,1)', ...
        %    'Colormap', gray)
        %quiver(rhou(:,:,1)', ...
        %    rhou(:,:,2)')
        %disp(Ue)
        colorbar('southoutside')
        grid off
        axis equal
        xlabel('i'); ylabel('j');
        drawnow
        % saveDat(Phi, sprintf('Phi%i', I))
    end
% ------- log 
end
disp("SIMULATION FINISHED")
r = 1;
end % ###################################################
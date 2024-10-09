% # net.m ###############################################
% Programm: 		    	simpleFlowSolver2D
% Content:					poisson solver
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% function that create kernel for solving the poisson eq.
% =======================================================

function A = poissonKernel2D(A, BCs) % ##################
% A::Matrix of size (Ima*Jma,Ima*Jma), initialized w 0s
% BCs::Vector of length 4, value of BC:
%   - 0 --> Dirichlet: specifying boundary value = 0
%   - 1 --> von Neumann: specifying boundary gradient = 0
%
%   _____ BC(4) _____
%  |                 |
%  |                 |
%  |                 |
% BC(3)             BC(1)
%  |                 |
%  |                 |
%  |_____       _____|
%         BC(2)
%
% ------- parameters ------------------------------------
global Ima Jma
%Ima = 4; Jma = 3;
BCKernel = zeros(Ima, Jma);
% ------- operation -------------------------------------
A = diag(-4*ones(1, Ima*Jma), 0);

BCKernel = BCs(1)*[zeros(Ima-1, Jma); ones(1, Jma)] +...
           BCs(2)*[ones(Ima, 1)  zeros(Ima, Jma-1)] +...
           BCs(3)*[ones(1, Jma); zeros(Ima-1, Jma)] +...
           BCs(4)*[zeros(Ima, Jma-1)  ones(Ima, 1)];
%heatmap(BCKernel)

for j = 1:Jma
    for i = 1:Ima
        ind = i + Ima*(j-1);
        A(ind,:) = funcNotZero( ...
            A(ind,:), @(x) x +BCKernel(i,j) );
    end
end

A = A + diag(ones(1, Ima*Jma-1), +1) + ...
        diag(ones(1, Ima*Jma-1), -1) + ...
        diag(ones(1, Ima*(Jma-1) ), +Ima) + ...
        diag(ones(1, Ima*(Jma-1) ), -Ima) ;
for k = 1:Jma-1
    A(Ima*k, Ima*k +1) = 0;
    A(Ima*k +1, Ima*k) = 0;
end

% for i = 1:Ima*Jma
%     for j = 1:Ima*Jma
%         A(i,j) = funcNotZero(A(i,j), @(x) 1/x);
%     end
% end

end % ###################################################


function val = funcNotZero(val, func) % #################
% ----- apply specified function to vector or scalar val
if val == 0
    return
end
len = length(val);
if len == 1
    val = func(val);
else
    for i = 1:len
        val(i) = funcNotZero(val(i), func);
    end
end
end % ###################################################

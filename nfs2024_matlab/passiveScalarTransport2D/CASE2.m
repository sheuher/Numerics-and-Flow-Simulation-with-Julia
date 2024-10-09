% # CASE2.m #############################################
% Programm: 		    	passiveScalarTranport2D
% Content:					define case, initiate fields
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% initiate fields with the defined case for case 2
% =======================================================

% #######################################################
% ------- parameters: transported field initialization --
Phi = zeros(ImaAll, JmaAll);
r   = int64(ImaAll/2);
for i = Ifi:Ila
    for j = Jfi:Jla
        if (i>r-r/7 && i<r+r/7 && j<r/5)
            Phi(i,j) = 1;
        end
    end
end
% ------- parameters: velocity field initialization -----
Delta = 1/Ima;
x = linspace( ((-ImaAll/2) +0.5) *Delta, ...
        ((ImaAll/2) -0.5) *Delta, ImaAll )';
y = linspace( ((-JmaAll/2) +0.5) *Delta, ...
        ((JmaAll/2) -0.5) *Delta, JmaAll )';
X = ones(length(x), 1) * x';
Y = y * ones(length(y), 1)';
%U = cos((X + 0.5)*2*pi) .* sin((Y + 0.5)*2*pi);
%V = - sin((X + 0.5)*2*pi) .* cos((Y + 0.5)*2*pi);
U = zeros(ImaAll, JmaAll);
V = zeros(ImaAll, JmaAll);
U(:,:) = 0;
V(:,:) = 3;
% ------- parameters: other fields initialization -------
dims =  [ImaAll, JmaAll];
Phie =  zeros(dims);
Phin =  zeros(dims);
Ue =    zeros(dims);
Vn =    zeros(dims);
FcX =   zeros(dims);
FcY =   zeros(dims);
FdX =   zeros(dims);
FdY =   zeros(dims);
PhiNew =zeros(dims);
% #######################################################
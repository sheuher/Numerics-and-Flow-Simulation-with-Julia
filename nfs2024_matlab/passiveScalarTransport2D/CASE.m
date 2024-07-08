% # CASE.m ##############################################
% Programm: 		    	passiveScalarTranport2D
% Content:					define case, initiate fields
% Author:                   S.H. Tey, June 2024
% Email:                    sheu.tey@stud.uni-due.de
% Matriculation Number:     3117568
% -------------------------------------------------------
% initiate fields with the defined case
% =======================================================

% #######################################################
% ------- parameters: transported field initialization --
Phi = zeros(ImaAll, JmaAll);
for i = Ifi:Ila
    for j = Jfi:Jla
        Phi(i,j) = sin( (j-nG-0.5)/Ima *10*pi );
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
U = cos((X + 0.5)*2*pi) .* sin((Y + 0.5)*2*pi);
V = - sin((X + 0.5)*2*pi) .* cos((Y + 0.5)*2*pi);
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
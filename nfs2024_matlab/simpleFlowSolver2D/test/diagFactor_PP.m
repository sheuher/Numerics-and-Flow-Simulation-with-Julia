function diagFactor_PP = diagFactor_PP(BCtype)
% returns positive diagonal factors of the pressure poisson equation 
%   BCtype(num_B=4) sets the B.C:
%       BCtype(1) --> Left --> j = Jfim
%       BCtype(2) --> Right --> j = Jlap
%       BCtype(3) --> Bottom --> i = Ifim
%       BCtype(4) --> Top --> i = Ilap
%   Types :
%       0 --> zero gradient
%       1 --> zero drichlet 
% dummy numbers exist in the diagFactor_pp in ghost cells

global ImaAll Ifi Ila 
global JmaAll Jfi Jla 

diagFactor_PP = ones(ImaAll,JmaAll)*4; 
modifier = zeros(4,1); % modifies the factors of each boundary;

for i= 1:4
    if (BCtype(i) == 0)
        modifier(i) = -1;
    elseif ( BCtype(i) == 1)
            modifier(i) = 0;
    end
end

% modify boundary neighbours
diagFactor_PP(:,Jfi) = diagFactor_PP(:,Jfi) + modifier(1); %Left
diagFactor_PP(:,Jla) = diagFactor_PP(:,Jla) + modifier(2); %Right
diagFactor_PP(Ifi,:) = diagFactor_PP(Ifi,:) + modifier(3); %Bottom
diagFactor_PP(Ila,:) = diagFactor_PP(Ila,:) + modifier(4); %Top

end


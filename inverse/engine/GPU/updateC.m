function C=updateC(C,ins,exc,dim,lr)
% Updates C
% Input:
%   C = M^-1 * K
%   exc = Vector of different points of excitation
%   dim = number of nodes
%   lr = tells us whether to look at left or right ventrical
% return:
%   updated C
step_delay=0;
if (step_delay >=0 && ins == step_delay+1);
    s = zeros(dim,1);
    if lr=='L'
        s(exc == 1 | exc==3)=0.5;
    elseif lr == 'R'
        s(exc == 2 | exc==4)=0.5;
    else
        s(exc>0.2)=0.5;
    end
 
    % Set wherever there has been excitation to 0 inside of C
    C(s == 0.5,:) = 0;  
end
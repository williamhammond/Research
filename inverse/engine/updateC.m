% Updates the gradient
% Input:
%   C = current gradient
%   ins = index 
%   exc = excitation
%   dim = number of nodes
% Return:
%   Updated C 
function C=updateC(C,ins,exc,dim)
step_delay=0;
if (step_delay >=0 && ins == step_delay+1);
    processC(C,exc,dim,'R');
end
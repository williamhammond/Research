% Gets Current Stimuli
% Input:
%   ins = index of stimulus array
%   dim = number of nodes
%   sti_l = stimulus array of left ventrical 
%   exc = excitation array
% Return:
%   New stimulus array
function sti = getCurrentSti(ins,dim,sti_l,exc)

sti=zeros(dim,1);

if sti_l(ins) == 1
    sti(exc == 1 | exc==3)=0.5;
end

% if sti_r(ins) == 1
%     sti(exc == 2 | exc==4)=0.5;
% end
%     

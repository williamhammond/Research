function[C, sti_l,exc]= setexcitation (path,C,dim,Ntime,TimeSeq)
%%in this code, the excitation only consider the normal excitiation,
% no-delay and no LBBB and RBBB
% If face LBBB or RBBB or other delay, the code need to add more options.
% This function updates the C matrix given the excitation array and builds
% and returns a stimulation arrray
% Input:
%   path = path to heart.exc
%   C = M^-1*k matrix
%   dim = number of nodes
%   Ntime = the number of steps in the time series
%   TimeSeq = time series ofin this code, the excitation only consider the normal excitiation, the propagation
% Return:
%   C = updated C matrix
%   sti_l = stimulus vector for left ventrical
%   exc = excitation array
sti_l=zeros(1,Ntime);

filename = strcat(path,'heart.exc');
fid = fopen(filename,'rb');
source = fread(fid,'int');
fclose(fid);

num_exc = source(1);
exc = source(2:end); 

duration = input('Please input the excitation duration (ms):','s');
duration_exc = str2double(duration);

step_delay = 0;
C=processC(C,exc,dim,'L');

% sti_l only needs value whenever stimulation is occuring
for i=1:Ntime
    if duration_exc >= TimeSeq(i)
        sti_l(i)=1;
    end
end
    


% if is_pacing == 0
%     %********determine initial conditions, by .exc & focus.seq*************
%    
% 
%     num_source = source(1);
%     source = source(2:end);   %important--to get rid of first number
%     %normal, while rv is a littelb bit delayed than lv
%     if delay>0
%         s(source == 1 | source == 3) = 0.5;    % lv first
%     end
%     %normal, simultaneous l/rv
%     if delay ==0
%         s(source == 1 | source == 3 | source ==2 | source == 4 ) = 0.5;
%     end
%     %rbbb
%     if delay == -1
%         s(source == 1 | source == 3) = 0.5;
%     end
%     %lbbb
%     if delay == -2
%         s(source ==2 | source == 4 ) = 0.5;
%     end
% 
% end
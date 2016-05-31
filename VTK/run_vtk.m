%% File
%   Runs VTK sequence generation code for Hackathon (KIT) datasets
%% Author
%   William Hammond
%% Date
%   4/4/2016
%% Implementation
path_in = '/home/wth4280/Research/TVR/Matlab/Output';

U = zeros(2223,264);
%% Put all signals into a sequence
for i = 1:264
    file = fullfile(path_in,strcat('outputx', num2str(i),'.bin'));
    
    fid = fopen(file,'rb');
    xi = fread(fid,[2223,inf],'double');
    fclose(fid);
    
    U(:,i) = xi;
end
write_vtk_sequence('v',U,1,'HackathonInverse');
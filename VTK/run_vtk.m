%% File
%   Runs VTK sequence generation code for Hackathon (KIT) datasets
%% Author
%   William Hammond
%% Date
%   4/4/2016
%% Implementation
load('Simulation_01_SEPTUMCENTER.mat');
U = Simulation_01_SEPTUMCENTER.potvals;
write_vtk_sequence('v',U,1,'output');
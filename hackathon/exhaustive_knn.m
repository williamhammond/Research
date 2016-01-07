% File: exhaustive_knn.m
% Author: William Hammond
% Date: 01/07/2016
% Description:  This code applies the k-nn algorithm on a set of synthetic 
%               TMP signals based off a single patient. The patient's 
%               geometry was used to exhaustively simulate the TMP 
%               signals treating every possible node in the heart as the 
%               source of the signal. 

tic;
% Path to directory holding all of the TMP signals 
dataPath = '../data/processed/Simulation/1898/Output/forward/';
% The number of nodes
dim = 1898;
% The number of timesteps we used 
nTime = 10153;

% Read in the ground truth data
i = 1;
fileName = sprintf('output%d.bin', i);
resultPath = strcat(dataPath,fileName);
fid = fopen(resultPath,'r');
target = fread(fid, [nTime,dim], 'double');
fclose(fid);

bestErr = inf;
bestNode = 0;
for i = 2:10
    % Read in output
    fileName=sprintf('output%d.bin', i);
    resultPath = strcat(dataPath,fileName);
    fid = fopen(resultPath, 'r');
    curr = fread(fid, [nTime,dim], 'double');
    fclose(fid);
    
    newErr = immse(target,curr);
    if newErr < bestErr
        bestErr = newErr;
        bestNode = i;
    end
end
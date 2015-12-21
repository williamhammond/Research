function [U_t,Ru,samV]=dynamicmodel(C,meanX,samV,Px,i,time,sti_l,exc)

%%The function is used to calculate the dynamic estimation
%   runs the Volumetric Myocardial TMP Activity Model
% Input
%   C = is equal to -M^-1 * K
%   meanX = 
%   samV = Sample of the current recovery
%   Px = Initial potential?
%   time = Time series of the propagation
%   time_inverse = Time series for the inverse
%   sti_l = Boolean array that tells us when to stimulate? 
%   exc = Array of points that are being stimulated?
%   noiCov = 

dim_x = length(meanX);
par = 0.15*ones(dim_x,1);
% Generate the weights
[samW, nu]=samWeiGen(dim_x);
% Generate the sample for TMP
samX=samGen(dim_x,meanX,Px,nu);
[samX,samV]=staSamPropar_gpu(C,samX,samV,par,i,time,exc,sti_l);

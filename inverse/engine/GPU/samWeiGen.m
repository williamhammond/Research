function [samW, nu]=samWeiGen(dim_x)
% This function generates a vector of weights. 
% Input:
%   dim_x = dimension of the current sample mean
ALPHA = 0.5;
BELTA = 2.0;
KAPA = 1.0;

lamda = (ALPHA*ALPHA*(dim_x+KAPA)-dim_x);
nu = sqrt(dim_x+lamda);

samW(1)=lamda/(dim_x+lamda);
samW(2)=1/(2*(dim_x+lamda));
samW(3)=samW(1)+1-ALPHA*ALPHA*BELTA;
samW(4)=samW(2);

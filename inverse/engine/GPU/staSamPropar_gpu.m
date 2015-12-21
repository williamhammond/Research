% Start Sample Propagation
% This function  
% Input:
%   C - M^-1 * K
%   samX - Sample of U (action potential)
%   samV - Sample of V (recovery current)
%   par - Parameter (a in our case 0.15)
%   k - The time index of the inverse
%   time - The time series for the forward model
%   time_inverse - The time series for the inverse
%   exc - Excitation array
%   sti_l - Stimulus array for the left ventrical

function [U,V]=staSamPropar_gpu(C,U,V,par,k,time,exc,sti_l)

    % Prepare data to be used on the GPU 
    [dim, dim_sam]=size(U);
    U=gpuArray(U);
    V=gpuArray(V);
    par=gpuArray(repmat(par,1,dim_sam));
    C=gpuArray(C);

    % Get current stimuli
    s=getCurrentSti(k,dim,sti_l,exc);
    % Update C 
    C=updateC(C,k,exc,dim,'L');

    ta = time(k);
    tb =  time(k+1);
    dt = tb - ta;

    [U,V]=samPro_RKA(U,V,par,s,dt,C);

    U=gather(U);
    V=gather(V);

    %fprintf('step = %d\n',count);
end
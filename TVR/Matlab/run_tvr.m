%%
% File: 
%   run_tvr.m   
%
% Author:
%   William Hammond 
%
% Description:
%   TV method to find the excitation location based off of KIT hackathon 
%   data. In order to run this, download Jinjja's TVR directory and
%   and execute this file on the same path as L2TV_region_loc.m
%
% Input:
%   path      - file path to the geometry, transfer matix, and mask
%   timeSteps - Number of time instances to generate
%
% Output:
%   Generates a matlab file for the potentials at each pacing site as well
%   as individual binary files for each gradient and potential at each time
%   instance
%%
function run_tvr(path, timeSteps)
    %% Data loading
    %%load data (include transfer matrix H, heart cor, input BSP and guide pacing location-carto map)

    %Transfer matrix
    file = strcat(path,'Trans.bin');
    fid=fopen(file,'rb');
    H=fread(fid,[2002,inf],'double');
    fclose(fid);

    % Mask
    file = strcat(path,'mask.idx');
    fid = fopen(file,'rb');
    mask = fread(fid,'int');
    fclose(fid);

    mask = logical(mask);
    H = H(mask,:);

    [n,m]=size(H);

    file = strcat(path,'deltasort.bin');
    fid=fopen(file,'rb');
    D=fread(fid,[3*m,m],'double');
    fclose(fid);

    %heart cor
    file = strcat(path,'heart.cor');
    fid = fopen(file,'rb');
    Cor = fread(fid,[3,inf],'double');
    fclose(fid);

    pacing_sites = ['SEPTUMCENTER'; 'RVPOSTERIOR '; 'RVANTERIOR  '; 'LVLATEPI    '; 'LVLATENDO   '; 'LVLAT       '; 'LVANTERIOR  ';'LVAPEX      '];
    pacing_sites = cellstr(pacing_sites);

    for i = 1:length(pacing_sites)
        site = pacing_sites(i);
        site = site{1};
        file = fullfile('Input', strcat('bsp_', site, '.mat'));

        % Evaluate as string to account for different bsp names
        BSP = load(file);
        exp = strcat('BSP.', strcat('bsp_', site));
        BSP = eval(exp);

        result = zeros(m,timeSteps); 

        for t = 1:timeSteps
            bsp = BSP(:,t);
            %% TV method 
            [x,dx] = TVR_method(path,bsp,H,D);

            result(:,t) = x;

            file_out = strcat('outputx_',site, num2str(t), '.bin');
            fid = fopen(strcat('Output/', file_out), 'wb');
            fwrite(fid,x,'double');
            fclose(fid);

            file_out = strcat('outputdx_', site, num2str(t), '.bin');
            fid = fopen(strcat('Output/', file_out), 'wb');
            fwrite(fid,x,'double');
            fclose(fid);

        end
        file_out = strcat('Output/', 'inverse_potval_', site, '.mat');
        save(file_out, 'result');
    end
end

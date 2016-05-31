%%
% File: 
%   L2TV_region_loc.m   
%
% Author:
%   Jingja Xu 
%
% Description:
%   TV method to find the excitation location  
%%
clear all
close all
clc
%% Data loading
%%load data (include transfer matrix H, heart cor, input BSP and guide pacing location-carto map)
path='/home/wth4280/Research/data/Hackathon/Simulation/2223/Input/';
timeSteps = 40;

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
    file_out = strcat('inverse_potval_', site, '.mat');
    save(file_out, 'result');
end

%% postprocessing: data-plot and find the excitation location
%figure, bar(x);  ylim([min(x) max(x)]); title('reconstructed signal x ');
% [val,idx]= max(x);
% Erro=norm(Cor(:,idx)-carto_pl(:,cartoid),2);
%thr =0.01;
%X_r=Cor(1,:);    i = 1
%Y_r=Cor(2,:);
%Z_r=Cor(3,:);
%figure;
%plot3(X_r,Y_r,Z_r,'g.');
%hold on
%lx=X_r(x>thr);
%ly=Y_r(x>thr);
%lz=Z_r(x>thr);
%plot3(lx,ly,lz,'b*');
%hold on
%ltx=carto_pl(1,cartoid);
%lty=carto_pl(2,cartoid);
%ltz=carto_pl(3,cartoid);
%plot3(ltx,lty,ltz,'r*');

% Save data



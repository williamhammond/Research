%test period ********************
%function state = staPro(par,dim,sti_amp, sti_per, sti_num,sti_start,sti_per2)
%test parameter influence*******
%function [ot,os] = staPro(s,dim,par)
%compute the state***********
%type = 0, using RK; type = 1, using ode
 function [t,state] = staPro(dim,a,ready,model,type,delay,is_pacing,varargin)
% % compute the state********************************************
% dim = dimension of state vector
%  a: =0.15 for 'P' / 0.13 for 'M'
% ready = 1 with Trans_state, else = 0
% model = 'P','O','M'
% type = 0 with RK / 1 with ode82
% delay = -1, RBBB / -2, LBBB/ 0 no delay / >0 dela
path_geo = '/mnt/hgfs/lxwast/Research/Data/HalifaxEP/M4_RX/Processed/Final/2230/';
%Cube/Final/1331/';Auckland/Auckland_f02_2_5/Processed/Results/Final/836/';PhysioNet/case4/Processed/Results/GenericGeometry/Final/1839/';

path_in = '/mnt/hgfs/lxwast/Research/Data/HalifaxEP/M4_RX/Processed/Simulation/2230/Input/N/';
%PhysioNet/case2/Processed/Results/Simulation/1373/Input/';Auckland/Auckland_f02_2_5/Processed/Results/Simulation/836/Input/Torso370/';
%Cube/Simulation/1331/Input/'4PhysioNet/case4/Processed/Results/GenericGeometry/Simulation/1839/Input/';
path_out = '/mnt/hgfs/lxwast/Research/Data/HalifaxEP/M4_RX/Processed/Simulation/2230/Output/N/';
%PhysioNet/case4/Processed/Results/Simulation/1373/Output/Test/';Auckland/A
%uckland_f02_2_5/Processed/Results/Simulation/836/Output/Torso370/';PhysioNet/case4/Processed/Results/GenericGeometry/Simulation/1839/Output/';
%ation/1331/Output/';
a_mid =  0.15; %0.1;%13;%15;
a_epi = 0.15;%0.17% 3;%7;
a_endo = 0.15;%0.14%3;%14;
a_isc = 0.3;
a_inf = 0.5;

if isempty(varargin) ~= 1
    aab = varargin{1};
end
if(ready==1)
    filename = strcat(path_in,'Trans_state.bin');
    fid = fopen(filename,'rb');
    s = fread(fid,[dim,inf],'double');
    fclose(fid);
else
    filename = strcat(path_in,'Mass.bin');
    fid = fopen(filename,'rb');
    M = fread(fid,[dim,inf],'double');
    fclose(fid);
    filename = strcat(path_in,'Stiff.bin');
    fid = fopen(filename,'rb');
    K = fread(fid,[dim,inf],'double');
    fclose(fid);
    s = -inv(M)*K;
    filename = strcat(path_in,'Trans_state.bin');
    fid = fopen(filename,'wb');
    fwrite(fid,s,'double');
    fclose(fid);
end

X0 = zeros(2*dim,1);
%fid = fopen('./input/auckland1720/auckland_source.txt');

if is_pacing == 0
    %********determine initial conditions, by .exc & focus.seq*************
    filename = strcat(path_geo,'heart.exc');
    fid = fopen(filename,'rb');
    source = fread(fid,'int');
    fclose(fid);

    num_source = source(1);
    source = source(2:end);   %important--to get rid of first number
    %normal, while rv is a littelb bit delayed than lv
    if delay>0
        X0(source == 1 | source == 3) = 0.5;    % lv first
    end
    %normal, simultaneous l/rv
    if delay ==0
        X0(source == 1 | source == 3 | source ==2 | source == 4 ) = 0.5;
    end
    %rbbb
    if delay == -1
        X0(source == 1 | source == 3) = 0.5;
    end
    %lbbb
    if delay == -2
        X0(source ==2 | source == 4 ) = 0.5;
    end

else
    filename = strcat(path_geo,'Pacing.map');
    fid = fopen(filename,'rb');
    foci = fread(fid,'double');
    fclose(fid);
    X0(find(foci~=0)) = 0.5;

end

%if foci(1)~=0
%    X0(foci(2:length(foci))) = 0.5;
%end
filename = strcat(path_in,'Init_u.bin');
fid = fopen(filename,'wb');
fwrite(fid,X0(1:dim),'double');
fclose(fid);

%**********determine parameters by ischemia.seq*************** 
par = a_mid*ones(dim,1);
% filename = strcat(path_geo,'heart.cell');
% fid = fopen(filename,'rb');
% cell= fread(fid,'double');
% fclose(fid);
% par(find(cell == 0)) = a_epi;
% par(find(cell == 1)) = a_endo;

% filename = strcat(path_geo,'heart_isc_normal.map');
% fid = fopen(filename,'rb');
% isc = fread(fid,'double');
% fclose(fid);
% if isempty(find(isc~=0)) ~= 1         % not empty
%     par(find(isc == 2)) = a_isc;
%     par(find(isc == 3)) = a_inf;
% end
% filename = strcat(path_in,'Init_par_normal.bin');
% fid = fopen(filename,'wb');
% fwrite(fid,par,'double');
% fclose(fid);

s(X0 == 0.5,:) = 0;

t = [];
sta = [];
%using RK, for the sake of testing the correctness of C++ code
if type == 0
    filename = strcat(path_in,'time');
    switch model
    case 'M'
        filename = strcat(filename,'_M');
    case 'P'
        filename = strcat(filename,'_P');
    otherwise
        filename = strcat(filename,'_O');
    end
    filename = strcat(filename,'.bin');
    fid = fopen(filename,'rb');
    t = fread(fid,'double');
    fclose(fid);
    step = length(t);
    
    state = zeros(dim*2,step);
    state(:,1) = X0;
    par = a*ones(dim,1);

    if(delay>0)
        tmp = find(t <= delay);
        step1 = length(tmp);
        for i = 2:step1
            delta_t = t(i) - t(i-1);
            [state(1:dim,i),state(dim+1:end,i)] = RK(state(1:dim,i-1),state(dim+1:end,i-1),delta_t,s,par);
        end
        X0 = state(1:dim,step1);
        if(delay>0)
            X0(source == 2 | source == 4) = 0.5;
            s(source==2 | source == 4,:) = 0;
        end
        state(1:dim,step1) = X0;
        for i = step1+1:step
            delta_t = t(i) - t(i-1);
            [state(1:dim,i),state(dim+1:end,i)] = RK(state(1:dim,i-1),state(dim+1:end,i-1),delta_t,s,par);
        end
    else
        for i = 2:step
            delta_t = t(i) - t(i-1);
            [state(1:dim,i),state(dim+1:end,i)] = RK(state(1:dim,i-1),state(dim+1:end,i-1),delta_t,s,par);
        end        
    end
        
else
    % use ode***********************************
    opt=odeset('NormControl','on');
    switch model
        case 'M'
            if (delay>0)
                [t1,state1] = ode45(@f_m,[0,delay],X0,opt,s,dim,par);
                X0 = state1(end,:);
                X0 = X0';
                X0(source==2| source == 4) = 0.5;
                s(source == 2| source == 4,:) = 0;
                [t2,state2] = ode45(@f_m,[delay,300],X0,opt,s,dim,par);
                state = vertcat(state1(1:length(t1)-1,:),state2);
                t = vertcat(t1(1:length(t1)-1),t2);
            else
                [t,state] = ode45(@f_m,[0,100],X0,opt,s,dim,par);
%                 filename = strcat(path_out,'TMP.bin');
%                 fid = fopen(filename,'a+b');
%                 for i=1:200/200
%                     display(i);
%                     [ttmp,stmp] = ode45(@f_m,[(i-1)*200,i*200],X0,opt,s,dim,par);
%                     X0 = stmp(end,:);
%                     t = vertcat(t,ttmp(1:end-1));
%                     fwrite(fid,stmp(1:end-1,1:dim)','double');
%                 end
%                 fclose(fid);
%                 state = stmp;

            end
                
        case 'P'
            [t,state] = ode45(@f_p,[0,175],X0,opt,s,dim,par);
%         filename = strcat(path_out,'TMP.bin');
%         fid = fopen(filename,'a+b');
%        % for i=1:250/250
%             display(i);
%             [ttmp,stmp] = ode45(@f_p,[(i-1)*250,i*250],X0,opt,s,dim,par);
%             X0 = stmp(end,:);
%             t = vertcat(t,ttmp(1:end-1));
%             fwrite(fid,stmp(1:end-1,1:dim)','double');
%         end
%         fclose(fid);
%         state = stmp;

    otherwise
            [t,state] = ode45(@f,[0,130],X0,opt,s,dim,par);
    end
%    state = state';
end

state = state(:,1:dim)';
if type == 0
filename = strcat(path_out,'TMP_RK');
else
filename = strcat(path_out,'TMP_mat');
end
switch model
case 'M'
    filename = strcat(filename,'_M');
case 'P'
    filename = strcat(filename,'_P');
otherwise
    filename = strcat(filename,'_O');
end
filename = strcat(filename,'.bin');
fid = fopen(filename,'wb');
fwrite(fid,state,'double');
fclose(fid);

filename = strcat(path_in,'time');
switch model
case 'M'
    filename = strcat(filename,'_M');
case 'P'
    filename = strcat(filename,'_P');
otherwise
    filename = strcat(filename,'_O');
end
filename = strcat(filename,'.bin');
fid = fopen(filename,'wb');
fwrite(fid,t,'double');
fclose(fid);
end


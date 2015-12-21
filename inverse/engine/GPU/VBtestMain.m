% This code may need to adjust parameters: noiCov, Eps, Alpha_s 
% input and output path
% path='../../data/PhysioNet/Case3/Processed/Simulation/Input/';
% path2='../../data/PhysioNet/Case3/Processed/Simulation/Output/';

path='/home/wth4280/Documents/inverse/data/processed/Simulation/1898/Input/';
path2='/home/wth4280/Documents/inverse/data/processed/Simulation/1898/Output/';

% transfer matrix H 
file = strcat(path,'Trans.bin');
fid=fopen(file,'rb');
H=fread(fid,[120,inf],'double');
fclose(fid);

[n,dim]=size(H); % get the size of input and output

%Transtate for dynamic model
filename = strcat(path,'Trans_state.bin');
fid = fopen(filename,'rb');
C_init = fread(fid,[dim,inf],'double');
fclose(fid);

% read measurement datat
% file = strcat(path2,'BSP.bin');
% fid=fopen(file,'rb');
% BSP=fread(fid,[n,inf],'double');
% fclose(fid);
% 
% % read truth -- TMP
% file = strcat(path2,'TMP.bin');
% fid=fopen(file,'rb');
% TMP=fread(fid,[n,inf],'double');
% fclose(fid);

file = strcat(path,'time_P.bin');
fid=fopen(file,'rb');
TimeSeq=fread(fid,'double');
fclose(fid);
Ntime = length(TimeSeq);


duration = input('Please input the excitation duration (ms):','s');
duration_exc = str2double(duration);

% Matrix holding all 
excMat = eye(dim);
%Other initial setting

fprintf('Beginning Parell Run...\n');
for i = 1:2
    exc = excMat(i,:);
    C=processC(C_init,exc,dim,'L');
        
    sti_l=zeros(1,Ntime);
    % sti_l only needs value whenever stimulation is occuring
    for j=1:Ntime
        if duration_exc >= TimeSeq(i)
            sti_l(j)=1;
        end
    end
    
    U = zeros(dim,1);
    V = zeros(dim,1);
    par = 0.15*ones(dim,1);
   
    QU = zeros(Ntime,dim);

    for j = 1:Ntime - 1 
        [U,V]= staSamPropar_gpu(C,U,V,par,j,TimeSeq,exc,sti_l);
        QU(j,:) = U;
    end
    
    fileName = sprintf('output%d.bin', i);
    resultPath = strcat(path2,'forward/',fileName);
    parsave(resultPath, QU);
end

for i = 1:Ntime - 1

end




% matlabpool close
fprintf('Finish\n');

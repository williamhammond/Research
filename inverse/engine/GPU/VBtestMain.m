% This code may need to adjust parameters: noiCov, Eps, Alpha_s 
% input and output path

path='../../../data/processed/Simulation/1898/Input/';
path2='../../../data/processed/Simulation/1898/Output/';

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

file = strcat(path,'time_P.bin');
fid=fopen(file,'rb');
TimeSeq=fread(fid,'double');
fclose(fid);
nTime = length(TimeSeq);

duration = input('Please input the excitation duration (ms):','s');
duration_exc = str2double(duration);

excMat = eye(dim);

fprintf('Beginning Parell Run...\n');
for i = 1:1
    exc = excMat(i,:);
    C=processC(C_init,exc,dim,'L');
        
    sti_l=zeros(1,nTime);
    % sti_l only needs value whenever stimulation is occuring
    for j=1:nTime
        if duration_exc >= TimeSeq(j)
            sti_l(j)=1;
        end
    end
    
    U = zeros(dim,1);
    V = zeros(dim,1);
    par = 0.15*ones(dim,1);
   
    QU = zeros(nTime,dim);

    for j = 1:nTime - 1 
        [U,V]= staSamPropar_gpu(C,U,V,par,j,TimeSeq,exc,sti_l);
        QU(j,:) = U;
    end
%     i = 1;
%     fileName = sprintf('output%d.bin', i);
%     resultPath = strcat(path2,'forward/',fileName);
%     parsave(resultPath, QU);
end
% matlabpool close

conn = database('/home/wth4280/Documents/Research/data/processed/SignalAW.db','','','org.sqlite.JDBC','URL');

fprintf('Finish\n');

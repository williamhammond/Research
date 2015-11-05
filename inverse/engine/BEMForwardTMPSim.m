clear
% setup the path where the files are located
% dont forget to setup matlab path to this path
path = '/home/will/Documents/Research/inverse/data/processed/Simulation/1898/Input/';
numSur= 120;
numVol= 1898;

% get K
fileK = strcat(path,'BEMK.bin'); 
fidK=fopen(fileK,'rb');
K=fread(fidK,[numSur,inf],'double');

% get Ka
fileKa = strcat(path,'Ka.bin'); 
fidKa=fopen(fileKa,'rb');
Ka=fread(fidKa,[1,inf],'double');

% get M
fileM = strcat(path,'MFree.bin'); 
fidM=fopen(fileM,'rb');
M=fread(fidM,[numVol,inf],'double');

% calcualte Trans.bin
% number of rows: numSur, number of columns : numVol
K = K.';
K = vertcat(K,Ka);
M = M.';
Trans = inv(ctranspose(K)*K)*ctranspose(K)*M;

%write trans to file

fileTrans = strcat(path,'Trans.bin');

fidTrans=fopen(fileTrans,'wb');
fwrite(fidTrans,Trans,'double');

fclose(fidK);
fclose(fidKa);
fclose(fidM);
fclose(fidTrans);
%*******************************
%calculate Trans_state.bin


%n = numVol;

%get M
fileMass = strcat(path,'Mass.bin'); 
fidMass=fopen(fileMass,'rb');
Mass=fread(fidMass,[numVol,inf],'double');

%get K
fileStiff = strcat(path,'Stiff.bin');
fidStiff=fopen(fileStiff,'rb');
KStiff=fread(fidStiff,[numVol,inf],'double');


% computer inv(M)*K
Trans_state=-pinv(Mass)*KStiff;

%erite to file
fileTransState = strcat(path,'Trans_state.bin');
fidTransState=fopen(fileTransState,'wb');
fwrite(fidTransState,Trans_state,'double');
fclose(fidTransState);
fclose(fidMass);
fclose(fidStiff);



% generate parameter.bin file

para = ones(numVol,1)*0.15;
fileP = strcat(path,'parameter.bin');

fidP=fopen(fileP,'wb');
fwrite(fidP,para,'double');

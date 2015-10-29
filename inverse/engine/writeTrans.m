function writeTrans
%path = '/Users/jd1336/Work/local_gccis-casci-cbl/Data/HalifaxEP/DC/Processed/Simulation/2144/Input/';



path=input('Enter the path : ','s');
id=input('Enter the number of nodes : ');
numSur=input('Number of leads : ');
numVol = id;
%numSur =120;


% path = '';
% numSur= 120;
% numVol= 1405;

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
%
file3 = strcat(path,'Trans.bin');

fid=fopen(file3,'wb');
fwrite(fid,Trans,'double');
fclose(fid);

% plot Trans
figure()
stem3(Trans)
title('Plot of Transfer Matrix "H" ');
xlabel('lead position')
ylabel('nodes in heart')
zlabel('value of H')
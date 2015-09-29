function writeParameter
path=input('Enter the path : ','s');
n=input('Enter the number of nodes : ');
%path = '/Users/jd1336/Work/local_gccis-casci-cbl/Data/HalifaxEP/AW/Processed/Simulation/1405/Input/';
numofnode=n;
para = ones(numofnode,1)*0.15;

file3 = strcat(path,'Parameter.bin');

fid=fopen(file3,'wb');
fwrite(fid,para,'double');
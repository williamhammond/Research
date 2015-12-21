function writeTransState
path=input('Enter the path : ','s');
n=input('Enter the number of nodes : ');

file1 = strcat(path,'Mass.bin'); 
fid=fopen(file1,'rb');
M=fread(fid,[n,inf],'double');
fclose(fid);

file2 = strcat(path,'Stiff.bin');
fid=fopen(file2,'rb');
K=fread(fid,[n,inf],'double');
fclose(fid);


Trans_s=-pinv(M)*K;


file3 = strcat(path,'Trans_state.bin');

fid=fopen(file3,'wb');
fwrite(fid,Trans_s,'double');

end
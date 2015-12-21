function outputfile3(path,i,tmp1,Btvx,time) 
   file1=strcat(path,num2str(i),'/tmp.bin');
   fid=fopen(file1,'wb');
   fwrite(fid,tmp1,'double');
   fclose(fid);
   
  
   file1=strcat(path,num2str(i),'/tmpBts.bin');
   fid=fopen(file1,'wb');
   fwrite(fid,Btvx,'double');
   fclose(fid);

   file1=strcat(path,num2str(i),'/time.bin');
   fid=fopen(file1,'wb');
   fwrite(fid,time,'double');
   fclose(fid);
  
%    file1=strcat(path,num2str(i),'/tmpdu.bin');
%    fid=fopen(file1,'wb');
%    fwrite(fid,du,'double');
%    fclose(fid);
%    
%    file1=strcat(path,num2str(i),'/tmpdxBtv.bin');
%    fid=fopen(file1,'wb');
%    fwrite(fid,dxBtv,'double');
%    fclose(fid);
%    
%    file1=strcat(path,num2str(i),'/tmpdxFtv.bin');
%    fid=fopen(file1,'wb');
%    fwrite(fid,dxFtv,'double');
%    fclose(fid);
   
  

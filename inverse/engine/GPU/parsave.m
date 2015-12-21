function parsave(fname, x)
    fid = fopen(fname,'w');    
    fwrite(fid,x,'double');
    fclose(fid);
end
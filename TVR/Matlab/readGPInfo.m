function readGPInfo(path,num_mfree)

%GP info needed to evaluate the L1 regularizer
% after read in the file, at each iteration
% the L1 norm is calculated as: Reg(gpts{i}.idx_dSha)
% =Reg(gpts{i}.idx_dSha) + 
% (gpts{i}.wei).*(gpts{i}.dTd).*norm(gpts{i}.dSha*U(gpts{i}.idx_sha))

% read in te weight for each gtps in integral


file = strcat(path,'gp.wei');
fid = fopen(file,'rb');
gpts_wei = fread(fid,'double');
fclose(fid);

gpts_num = length(gpts_wei);

gpts = cell(gpts_num,1); 

% read in the number and index of neighbors for each gtps
file = strcat(path,'gp_inf.num');
fid = fopen(file,'rb');
gpts_infnum = fread(fid,'int');
fclose(fid);


file = strcat(path,'gp_inf.seq');
fid = fopen(file,'rb');
gpts_inf = fread(fid,'int');
fclose(fid);

% read in the shape functions for each gtps
file = strcat(path,'gp.shape');
fid = fopen(file,'rb');
gpts_sha = fread(fid,'double');
fclose(fid);

idx_inf = 1;
idx_sha = 1;

for i=1:gpts_num
    inf = gpts_inf(idx_inf:idx_inf+gpts_infnum(i)-1);
    idx_inf = idx_inf + gpts_infnum(i);
    
    shape = gpts_sha(idx_sha:idx_sha+4*gpts_infnum(i)-1);
    idx_sha = idx_sha + 4*gpts_infnum(i);
    shape = reshape(shape,4,gpts_infnum(i));
    
    dSha = shape(2:4,:);

    dTd = dSha'*dSha;
    idx_dTd = ones(length(inf),1)*((inf'-1).*num_mfree) + inf*ones(1,length(inf));
    
    gpts{i} = struct('wei',{gpts_wei(i)},'dSha',{dSha},'idx_sha',{inf},'dTd',{dTd},'idx_dTd',{idx_dTd});
end


file = strcat(path,'gpts2.mat');

save(file,'gpts');
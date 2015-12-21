function  v = Calculate_v(cov_x,x,id)
% file=strcat(path,filename);
% load(file);
%sha= zeros(m,1);
global gpts;

gpts_num = length(gpts);
v = zeros(gpts_num,1);
%v2 = zeros(gpts_num,1);

parfor i=1:gpts_num
    U = x(gpts{i}.idx_sha);
    sha = gpts{i}.dSha;
    %sha2 = gpts{i}.dSha;
    %v(i) = trace(sha'*sha*(covx+U*U'));
    %v(i) = trace(sha2'*sha2*covx)+U'*sha'*sha*U
    if id == 1;
    v(i) =U'*sha'*sha*U;
    else
    covx = cov_x(gpts{i}.idx_sha,gpts{i}.idx_sha);
    v(i)=trace(sha'*sha*(covx+U*U'));%trace(sha'*sha*covx)+U'*sha'*sha*U;
    end
end 
    
end 
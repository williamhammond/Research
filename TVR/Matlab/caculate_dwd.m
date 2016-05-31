function [Regre]=caculate_dwd(path,U,m)


file=strcat(path,'gpts2.mat');
load(file);

Reg=zeros(m*m,1);
gpts_num = length(gpts);
for i=1:gpts_num
    dgpts=norm(gpts{i}.dSha*U(gpts{i}.idx_sha));
      if dgpts<1e-6
          dgpts=1e-6;
      end
   gpts1=(gpts{i}.wei).*(gpts{i}.dTd)./dgpts;
    id_dsha=gpts{i}.idx_dTd;
    id_dsha=id_dsha(:);
    gpts12=gpts1(:);
    Reg(id_dsha)= Reg(id_dsha)+ gpts12;
    
end 
Regre=reshape(Reg,m,m);

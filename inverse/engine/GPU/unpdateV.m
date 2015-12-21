function [samV,meanV] = unpdateV(samV,samW)

[dim_x,dim_sam]=size(samV);
meanV = samV(:,1)*samW(1);
meanV = meanV+sum(samV(:,2:end)*samW(2),2);%sum colum

for i = 1:dim_sam
    samV(:,i)=meanV;
end
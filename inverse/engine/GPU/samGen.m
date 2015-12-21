function samX=samGen(dim,meanX,Px,nu)
% This function generates a sample
% Input:
%   dim - the number of nodes
%   meanX - the sample mean
%   Px - the initial propagation? 
samX=zeros(dim,2*dim+1);
samX(:,1)=meanX;
tt = repmat(meanX,1,dim);
samX(:,2:dim+1)=tt+Px*nu;
samX(:,dim+2:end) =tt-Px*nu;
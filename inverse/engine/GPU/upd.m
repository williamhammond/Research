function [mean,P] =upd(sam, noiCov,samW)
[~,dim_sam]=size(sam);
mean = sam(:,1)*samW(1);
mean = mean+sum(sam(:,2:end)*samW(2),2);

sam = sam-repmat(mean,1,dim_sam);

P=sam(:,1)*sam(:,1)'*samW(3);
P=P+sam(:,2:end)*sam(:,2:end)'*samW(4);


P=P+noiCov;


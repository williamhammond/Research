% function Beta =  ExpBeta(r,mb,y,H,cov_x,x,n)
% Inv_Beta = r*(1/mb)+(1-r)/n/2*(norm(y-H*x)^2+trace(cov_x*H'*H));
% Beta = 1/Inv_Beta;
% end 

function [Beta,a,b] =  ExpBeta(a0,b0,y,H,cov_x,x,n)
% Inv_Beta = r*(1/mb)+(1-r)/n/2*(norm(y-H*x)^2+trace(cov_x*H'*H));
% Beta = 1/Inv_Beta;
a = a0+ n/2;
if b0==0
ivb = 1/2 *norm(y-H*x)^2 +1/2*(trace(cov_x*H'*H));
else
ivb = 1/b0+1/2 *norm(y-H*x)^2 +1/2*(trace(cov_x*H'*H));
end
b = 1/ivb;
Beta = a*(1/ivb); 


 end 
%  


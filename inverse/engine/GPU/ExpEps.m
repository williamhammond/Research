function [Eps,a,b] =  ExpEps(a0,b0,U,Ut,Us,Cov_Ut,Cov_U,Q,m)
% Inv_Beta = r*(1/mb)+(1-r)/n/2*(norm(y-H*x)^2+trace(cov_x*H'*H));
% Beta = 1/Inv_Beta;
a = a0+ m/2;
if b0==0
ivb = 1/2*(trace(U'*U)-2*U'*diag(Us)*Ut+Ut'*Q*Ut+trace(Cov_U)+trace(Cov_Ut*Q));
else
ivb = 1/b0+1/2*(trace(U'*U)-2*U'*diag(Us)*Ut+Ut'*Q*Ut+trace(Cov_U)+trace(Cov_Ut*Q));
end
b = 1/ivb;
Eps = a*(1/ivb); 


 end 
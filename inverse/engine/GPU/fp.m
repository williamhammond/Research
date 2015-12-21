function [u_tmp,v_tmp] = fp(u,v,s,Trans,par)



 k = 8;
 e = 0.01;

 
u_tmp = Trans*u + k*u.*(1-u).*(u-par)-u.*v+s;
v_tmp = -e*(k*u.*(u-par - 1) + v);
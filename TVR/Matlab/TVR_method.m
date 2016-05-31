function [x,dx] = TVR_method(path,bsp,A,D)
%Using TV- reweight method
%Input: path: data input path for gaussion point information.
%       bsp: input measurment information.
%       A: forward transfer matrix
%       D: Gradient matrix 
disp ('TV method')
IniX= TestX0_region(path,bsp,A);
dxIni = CoGradient_s(D,IniX,1);
m = length(IniX);
%%L2TV
x=IniX; 
dx1=dxIni;
   for sq=1:9
     Dwd=caculate_dwd(path,x,m);
     lamda_u=norm(A'*A,inf)/norm(Dwd,inf);
     Ax=A'*A+lamda_u*Dwd;
     bx=A'*bsp;
     U_k=pcg(Ax,bx,1e-7,50000);
     x=U_k;
     %R=[R x];
     dx= CoGradient_s(D,x,0);
     %norm(dx1-dx,2)
     if norm(dx1-dx,2)<0.001
       break;
     end
     dx1=dx;
   end 
dx= CoGradient_s(D,x,1);
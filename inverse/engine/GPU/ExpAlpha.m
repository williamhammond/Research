% function Alpha = ExpAlpha(r,ma,v,m,id)
% % r = a1/(a1+id*m);
% % ma = a1/b1;
% Inv_Alpha = r*(1/ma)+(1-r)*sum(sqrt(v))./(id*m);
% Alpha = 1/Inv_Alpha;
% end 


 function [Alpha,a,b] = ExpAlpha(a0,b0,v,m,id)
% % r = a1/(a1+id*m);
% % ma = a1/b1;
% Inv_Alpha = r*(1/ma)+(1-r)*sum(sqrt(v))./(id*m);
% Alpha = 1/Inv_Alpha;

 a = a0+id*m;
 if b0==0
      ivb =sum(sqrt(v));
 else
     ivb =1/b0+sum(sqrt(v));
 end
 b=1/ivb;
 Alpha = a*b;
 end 




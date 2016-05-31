function deltax = CoGradient_s(D,x,ig)

m=size(x,1);
Dx=D(1:m,:);
Dy=D(m+1:2*m,:);
Dz=D(2*m+1:end,:);

deltax=sqrt((Dx*x).^2+(Dy*x).^2+(Dz*x).^2);

%if(ig==1)
%figure
%subplot(2,1,2); bar(deltax);  ylim([min(deltax) max(deltax)]); title('delta X');
%end
% figure
% subplot(3,1,1); bar(Dx*tmp1);
% subplot(3,1,2); bar(Dy*tmp1);
% subplot(3,1,3); bar(Dz*tmp1);
% 
% figure
% subplot(3,1,1); bar(Dx*x);
% subplot(3,1,2); bar(Dy*x);
% subplot(3,1,3); bar(Dz*x);
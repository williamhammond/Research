
function x12=TestX0_region(path,bsp2,H)
%weight lap method caculate Xini
file = strcat(path,'heart.cor');
fid=fopen(file,'rb');
Cor=fread(fid,[3,inf],'double');
fclose(fid);


disp ('second order')
[W,B]=LSMN_mesh_gen(H,Cor,5); 
L=B*W;    
[U,sm,X] = cgsvd(H,L);
% figure
% lambda=l_curve(U,sm,bsp2);
lambda=0.5;
x12= tikhonov(U,sm,X,bsp2,lambda);

%figure
%bar(x12);  ylim([min(x12) max(x12)]); title('x second in noise');




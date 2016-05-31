function [ W B1 ] = LSMN_mesh_gen( H,cord,N )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% caculate W  
l=size(H,2);
w=zeros(l,1);
for i=1:l
w(i)=sqrt(H(:,i)'*H(:,i));
end

W=diag(w);

%% caculate B ( domain < 10 )

B= ones(l,l);
diff=zeros(1,l);
num=zeros(1,l);
for i=1:l
    for j=1:l
diff(j)=norm(cord(:,j)-cord(:,i));
    end
marsk = (diff<N);
num(i)=size(find(marsk==1),2);
B(i,:)=B(i,:).*marsk; 
end 
num2=diag(num+1);
% num2=diag(ones(1,l));
B1=(B-num2);



end


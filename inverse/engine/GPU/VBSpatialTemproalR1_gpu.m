function [TmpX,Ute,R] = VBSpatialTemproalR1_gpu(y,H,Ut_pre,R,iter)
%VB method for spatial and temporal decompostion 
%U = diag(Us)*Ut

%% Inital 

[n,m] = size(H);
I = ones(m,1);
%U = gpuArray(zeros(m,iter));
Us = rand(m,1);
Ut = Ut_pre;

%for mesurment noise
Eps = 1000;   %for sptial-temporal decomposition (U=diag(Us)*Ut)
ae0 = m;
be0 = Eps/ae0;

am0 = 2*n;   %for bsp = H*U
bm0 = n/trace(y'*y);
Beta =am0*bm0;

Cov_Ut=gpuArray(R);%diag(0*I);
v = Calculate_v(diag(0*I),Us,1);  %% the inital v for VB
%read gpu matrix

H=gpuArray(H);
IR = gpuArray(pinv(R));
% Us = gpuArray(Us);
% Ut = gpuArray(Ut);

%% VB method step 
%each step of VB iteration method update 
for k = 1: iter
%calculate U

    TempU = Beta* H'*H + Eps*eye(m);
    Cov_U = diag(I)/TempU;
%     [u,s,vd]=svd(TempU);
%     Ts =diag(s);
%     Ns = length(find(Ts>1e-9));
%     Cov_U = vd(:,1:Ns)*pinv(s(1:Ns,1:Ns))*u(:,1:Ns)';
    Y_U = Beta*H'*y+Eps*Us.*Ut;
    U =TempU\ Y_U;
 
%calculate Us
   %calculate Us within U = diag(Ut)*Us
 %Alpha_s initial set 1 or 10
    Alpha_s = ExpAlpha(m,10/m,v,m,1);
   
    Wv = dWd(m,gather(v));
    Wv = gpuArray(Wv);
    P = diag(diag(Ut*Ut'+Cov_Ut));
    TempUs = Eps.*P + Alpha_s.*Wv;
    Cov_Us = diag(I)/TempUs;
    Y_Us = diag(Eps*U*Ut');
    Us = TempUs\Y_Us;
%     Us(:,k+1) = pcg(TempUs,Y_Us,1e-6,50000);   
   
%calculate Ut
     Q = diag(Us)'*diag(Us) + diag(diag(Cov_Us)); 
     TempUt = Eps.*Q + IR;
     Cov_Ut = diag(I)/TempUt;
     Y_ut2= Eps.*diag(Us)'*(U-diag(Us)*Ut_pre);
     teX =TempUt\Y_ut2;
     Ut=teX+Ut_pre;
%      K=Eps.*diag(Us)'*Cov_Ut;
%      IR =(eye(m)-K*diag(Us))*IR;
%      yt=Eps*diag(Us)*U+(IR)*Ut_pre;
%      Ut = TempUt\yt;
    %Ut2= Ut_pre+Eps.*Cov_Ut*diag(Us(:,k+1))*(U(:,k+1)-diag(Us(:,k+1))*Ut_pre);
    
%Hyperparameters update
     Beta = ExpBeta(am0,bm0,y,H,Cov_U,U,n);
     Eps = ExpEps(ae0,be0,U,Ut,Us,Cov_Ut,Cov_U,Q,m);
     %Alpha_t = ExpAlpht(at0,bt0,Ut,Ut_pre,Cov_Ut,m);
     v = Calculate_v(gather(Cov_Us),gather(Us),0);

% 

end 

      
      
%% Output         
 Ute = gather(Ut); 
 TmpX = gather(2*U);
 R = gather(Cov_Ut);
      
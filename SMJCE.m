function [NMSE_SMJCE,NMSE_SMMV]=SMJCE(Yk,YkHermite,V,G,System,N_bpilot) % Algorithm 2  
 
NMSE_SMMV=0;
NMSE_SMJCE=0;
NM=System.NM;  
NL=System.NL;  
K=System.K;   
SNR=System.SNR;
Gr=System.Gr; 
N_npilot=System.N_npilot;
Dr=System.Dr;  

%%  Channelsubpace Estimation and Projection========================================
[~,Subpace]=CSMUSIC_AOD_Estimation( N_bpilot, Yk, YkHermite,System); % Subspace Estimation eq(26)
[~,length_ind]=size(Subpace);
YkProjected=zeros(N_bpilot,length_ind,K);
for ite=1:1:K
    YkProjected(:,:,ite)= YkHermite(:,:,ite)*Subpace*(Subpace'*Subpace)^(-1); %eq(28)
end 
noisecovariance= N_bpilot*trace((Subpace'*Subpace)^(-1))/(SNR)/(N_npilot);
G_hatall= zeros(NL,NM,K);
for itekkk=1:1:K
    Inedxk=CS_MMV_SOMP(YkProjected(:,:,itekkk),V'*Dr,noisecovariance);
    leftinverse=((V'*Dr(:,Inedxk))'*V'*Dr(:,Inedxk))^(-1)*(V'*Dr(:,Inedxk))';
    Xk=leftinverse*YkProjected(:,:,itekkk);
    if itekkk==1 
        X_curini= zeros(Gr,length_ind);
        X_curini(Inedxk,:)=Xk;
    end
    G_hat=Dr(:,Inedxk)*Xk*Subpace';
    G_hatall(:,:,itekkk)=G_hat;
    NMSE_SMMV=NMSE_SMMV+trace((G_hat-G(:,:,itekkk))'*(G_hat-G(:,:,itekkk)))/trace((G(:,:,itekkk))'*(G(:,:,itekkk)))/K; 
end 
%% MJCE
 [Alphaalll,Xall]=SMJCE_iteration(N_bpilot,V,G_hatall,X_curini,System,YkProjected,Subpace);
for itek=1:1:K
    G_hat=(Alphaalll((itek-1)*NL+1:(itek*NL),:))*Dr*Xall*Subpace';
    NMSE_SMJCE=NMSE_SMJCE+trace((G_hat-G(:,:,itek))'*(G_hat-G(:,:,itek)))/trace((G(:,:,itek))'*(G(:,:,itek)))/K;
end 
end
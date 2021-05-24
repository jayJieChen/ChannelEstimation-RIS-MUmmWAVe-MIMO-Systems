function [R_reduce,ChannelSubpace]=CSMUSIC_AOD_Estimation(N_bpilot,Yk,YkHermite, System)
 
z_bktilde=YkHermite;
z_bk=Yk;


NM=System.NM; % The number of antennas at the BS 
K=System.K;  % The numer of users
  

CovarianceZ=zeros(NM,NM);
for ite=1:1:K
    CovarianceZ=CovarianceZ+z_bk(:,:,ite)*z_bktilde(:,:,ite)./K;
end
[S_CovarianceZ,CovarianceZ_eigvalue,~]=eig(CovarianceZ);
[CovarianceZ_eigvalue_Trans, idexmaxtomin]=dsort(diag(CovarianceZ_eigvalue));
S_CovarianceZ_Trans=S_CovarianceZ;
for itej=1:1:length((idexmaxtomin))
    S_CovarianceZ_Trans(:,itej)=S_CovarianceZ(:,idexmaxtomin(itej));
end
CovarianceZ_eigvalue_Trans=abs(CovarianceZ_eigvalue_Trans);
%method 2
%p=Nb;
%N=observations:
MDLp=NM;
MDLN=K*N_bpilot;
MDL_value_cur=inf;
for MDLk=1:1:NM
    MDL_value_bef=MDL_value_cur;
    Term1=prod(CovarianceZ_eigvalue_Trans(MDLk+1:MDLp).^(1/(MDLp-MDLk)));
    Term2=sum(CovarianceZ_eigvalue_Trans(MDLk+1:MDLp))./((MDLp-MDLk));
    MDL_value_cur=-((MDLp-MDLk)*MDLN)*log((Term1/Term2))+0.5*MDLk*(2*MDLp-MDLk)*log(MDLN);
    if MDL_value_cur>MDL_value_bef
        break
    end
end


r=min(min(MDLk-1,NM),N_bpilot);%===!!===
ChannelSubpace=S_CovarianceZ_Trans(:,1:r);
R_reduce=sum(CovarianceZ_eigvalue_Trans(r+1:end)); 
end
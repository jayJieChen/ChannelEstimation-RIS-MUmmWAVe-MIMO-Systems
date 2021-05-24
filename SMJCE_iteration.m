function   [alphaAll,Xall]=SMJCE_iteration(N_bpilot,V,G_hatall,X_curini,System,YkProjected,Subpace);
NM=System.NM; % The number of antennas at the BS
NL=System.NL;  % The numer of reflective elements at the RIS
K=System.K;  % The numer of users
SNR=System.SNR;
Gr=System.Gr; 
N_npilot=System.N_npilot;
Dr=System.Dr;  

[~,length_ind]=size(Subpace);
lambdaD=0.1;
noisecovariance2=1/(SNR(1))/(N_npilot);
itemaxgra_Gradientsearch=4;
threshold_eps=1e-9;
%% First Estimated/ Joint Optimized 
Hestimated_cur=ones(NL,K);
for itek=1:1:K
    for itel=1:1:NL
        Hestimated_cur(itel,itek) =sum(G_hatall(itel,:,itek)./G_hatall(itel,:,1))./NM;
    end
end
Z_all_her=zeros(K*N_bpilot,length_ind);
V_all_her=zeros(K*N_bpilot,K*NL);
A_all=zeros(K*NL,NL);
for ite=1:1:K
    Z_all_her((ite-1)*N_bpilot+1:ite*N_bpilot,:)= YkProjected(:,:,ite);
    V_all_her((ite-1)*N_bpilot+1:ite*N_bpilot,(ite-1)*NL+1:ite*NL)=V';
    A_all((ite-1)*NL+1:ite*NL,:)= diag(Hestimated_cur(:,ite));
end
Z_all=Z_all_her';
V_all=V_all_her'; 
A_all_ini=A_all; 
lambda_ini=lambdaD/(log(Gr)*noisecovariance2);
lambda_all=diag(ones(1,K*N_bpilot));
for itee=1:1:K
    lambda_all((itee-1)*N_bpilot+1:itee*N_bpilot,(itee-1)*N_bpilot+1:itee*N_bpilot)=diag(lambda_ini.*ones(1,N_bpilot));  %/trace(noise(:,:,itee)'*noise(:,:,itee));%((trace((z_bktilde(:,:,itee)-V'*diag(Hestimated(:,itee)')*Dr*X_cur)'*(z_bktilde(:,:,itee)-V'*diag(Hestimated(:,itee)')*Dr*X_cur))));
end 
Xall= ones(Gr,length_ind);
obj_curout=sum(log(real(diag(Xall*Xall'))+threshold_eps))+trace((Z_all'-V_all'*A_all*Dr*Xall)'*lambda_all*(Z_all'-V_all'*A_all*Dr*Xall));

ind=1:1:NL;
ind=ind+NL*(ind-1);
for ite_sca=1:1:500
    obj_befout=obj_curout ;
    Abef_all=A_all;
    %% Compite X_cur
    obj_curinner=obj_befout;
    conver1=1e-6;
    for ite_scag=1:1:10000
        obj_befinner=obj_curinner;
        DiagnolaValue=diag(1./(real(diag(Xall*Xall'))+threshold_eps));
        Midterm1=(DiagnolaValue+(Dr'*A_all'*V_all*lambda_all*V_all'*A_all*Dr))^(-1);
        Xall=Midterm1*(Dr'*A_all'*V_all*lambda_all*Z_all');
        DiagnolaValue=diag(1./(real(diag(Xall*Xall'))+threshold_eps));
        Xall=(DiagnolaValue+Dr'*A_all'*V_all*lambda_all*V_all'*A_all*Dr)^(-1)*(Dr'*A_all'*V_all*lambda_all*Z_all');
        obj_curinner=real(sum(log(real(diag(Xall*Xall'))+threshold_eps))+trace((Z_all'-V_all'*A_all*Dr*Xall)'*lambda_all*(Z_all'-V_all'*A_all*Dr*Xall)));
        if (obj_befinner-obj_curinner)<conver1*1e-1
            break
        end
    end
    obj_curout=obj_curinner;
    
    %% Computde Derivate of A_ll
    for itegradient=1:1:itemaxgra_Gradientsearch
        if K==1
            flag=1;
            break;
        end
        
        for iteg=2:1:K
            vecZk=vec(YkProjected(:,:,iteg));
            Dmatrix11=kron((Dr*Xall).',V');
            Dmatrix11_ind=Dmatrix11(:,ind);
            flag=0;
            if rank(Dmatrix11_ind,1e-8)<NL
                flag=1;
                for itenu=1:1:NL
                    Abef_all((iteg-1)*NL+itenu,itenu)=A_all((iteg-1)*NL+itenu,itenu);%-DerivateA((iteg-1)*Nu+itenu,itenu)*(1e-4/(itegradient)/ite_sca^2)*abs(A_all_ini((iteg-1)*Nu+itenu,itenu)) ;%(1e-7/(itegradient)/ite_sca^2)*DerivateA((iteg-1)*Nu+itenu,itenu);
                end
            else
                Abef_all((iteg-1)*NL+1: iteg*NL,:)=diag( (Dmatrix11_ind'*Dmatrix11_ind)^(-1)*Dmatrix11_ind'*vecZk);
            end
            if  flag==1
                break
            end
        end
        obj_cur0=sum(log(real(diag(Xall*Xall'))+threshold_eps))+trace((Z_all'-V_all'*Abef_all*Dr*Xall)'*lambda_all*(Z_all'-V_all'*Abef_all*Dr*Xall));
        if obj_curout-obj_cur0>=0%real(trace(X_cur*X_cur'*DiagnolaValue))>real(trace(CovarianceZ_est'*((Dt_opt*DiagnolaValue_inv*Dt_opt')^(-1))*CovarianceZ_est))
            break
        end
        if  flag==1
            break
        end
    end
    
    if  flag==1
        break
    end
    
    if itegradient<itemaxgra_Gradientsearch
        A_all=Abef_all;
    end
    
    if abs(obj_curout-obj_cur0)<conver1
        break
    end
    obj_curout=obj_cur0;
end
alphaAll=A_all;
G_hat1=0;
G_hatini=0;
for itek=1:1:K
    G_hat1=sum(sum(abs((alphaAll((itek-1)*NL+1:(itek*NL),:))*Dr*Xall*Subpace').^2))+G_hat1;
    G_hatini=G_hatini+sum(sum(abs((A_all_ini((itek-1)*NL+1:(itek*NL),:))*Dr*X_curini*Subpace').^2));
end

if G_hat1>1.7*G_hatini %heuristically
    alphaAll=A_all_ini;
    Xall=X_curini;
end
end





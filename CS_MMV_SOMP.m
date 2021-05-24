function  Inedxk=CS_MMV_SOMP(z_bk2,Dmatrixo,yita)

[length01,length2]=size(Dmatrixo);
 
Rk=z_bk2;
AMplitude=zeros(length2,1);

Dmatrix=Dmatrixo;
for ii=1:1:length2
    AMplitude(ii)=sqrt(real(Dmatrixo(:,ii)'*Dmatrixo(:,ii)));%            norm(D_matxo(:,ii),2);
    if AMplitude(ii)<=1e-12
        AMplitude(ii)=1;
    else
        Dmatrix(:,ii)=Dmatrixo(:,ii)/AMplitude(ii);
    end
end

IdentiMatrix=diag(ones(length01,1));
InedxOmega_k=zeros(1,length2);
for iteSup=1:1:length2
    erv=sum(abs(Dmatrix'*Rk).^2,2);
    [~,jmax]=max(erv);
    InedxOmega_k(1,jmax(1))=1;
    Inedxk=find(InedxOmega_k(1,:)==1);
    it=0;
    while rank(Dmatrix(:,Inedxk))<length(Inedxk)&&it<20 %yita/;%yita%norm(vecRk,2)^2
        InedxOmega_k(1,jmax(1))=0;
        erv(erv==jmax)=0; 
        [~,jmax]=max(erv);
        InedxOmega_k(1,jmax(1))=1;
        Inedxk=find(InedxOmega_k(1,:)==1);
        it=it+1;
    end
    P_omega=Dmatrix(:,Inedxk)*(Dmatrix(:,Inedxk)'*Dmatrix(:,Inedxk))^(-1)*Dmatrix(:,Inedxk)'; 
    Rk=(IdentiMatrix-P_omega)*z_bk2;
    if real(trace(Rk'*Rk))<=yita %yita/;%yita%norm(vecRk,2)^2
        break;
    end
end

end

function LISpilotO=RIS_SequenceOptimization(System,N_bpilot)
Gr=System.Gr;
Dr=System.Dr;
NL=System.NL; 
xbins1=0:0.01:0.09;
LISpilot=sqrt(0.5)*(normrnd(0,1,NL,N_bpilot) + 1j*normrnd(0,1,NL,N_bpilot));%1-ny-nsignalzeros(n_signal,1);%
while rank(LISpilot)<N_bpilot
    LISpilot=sqrt(0.5)*(normrnd(0,1,NL,N_bpilot) + 1j*normrnd(0,1,NL,N_bpilot));%1-ny-nsignalzeros(n_signal,1);%
end
LISpilot=LISpilot./abs(LISpilot);
 

LISpilotO=LISpilot;
[CA1,CA0]=size(LISpilot);
A=Dr*Dr';
[Sx,Vx,Dx]=eig(A); %
[valuemaxtoimin, idexmaxtomin]=dsort(diag(Vx));
Sxtrans=Sx;
for itej=1:1:length((idexmaxtomin))
    Sxtrans(:,itej)=Sx(:,idexmaxtomin(itej));
end
CA2=length(valuemaxtoimin(valuemaxtoimin>0));
cur=inf;
for ite=1:1:1
    cuf=cur;
    for jj=1:1:CA0
        Tamma=Sxtrans'*LISpilot;
        VV=diag(valuemaxtoimin)*Tamma;
        Ej=CA0.*diag(valuemaxtoimin)-VV*VV'+VV(:,jj)*VV(:,jj)';
        [Sj,Vj,~]=eig(Ej); %
        [c1,~]=find(real(Vj)==(max(max(real(Vj)))));
        temvec=LISpilot(:,jj);
        for ii=1:1:CA2
            temvec(ii)=Sj(ii,c1(1))*sqrt(Vj(c1(1),c1(1)))/valuemaxtoimin(ii);
        end
        LISpilot(:,jj)=Sxtrans*temvec;
        LISpilot(:,jj)= LISpilot(:,jj)./abs(LISpilot(:,jj));
    end
    [~,cur,~]=Muturalance_correlance(Gr,LISpilot,Dr,xbins1); 
    if cur<cuf
        LISpilotO=LISpilot;
    else
        break
    end
end 
LISpilotO=LISpilotO./abs(LISpilotO); 
end

function  [F,H,GenieaidedAoA,GenieaidedAoD]=ChannelGenralize(System) ; 
 
NL=System.NL;
NM=System.NM;
Nf=System.Nf;
Nhk=System.Nhk;
K=System.K; 
AoAF=zeros(NL,Nf) ;
GenieaidedAoD=zeros(NM,Nf);
AoDHk=zeros(NL,Nhk,K);
GenieaidedAoA=zeros(NL,Nhk*Nf,K); 
% BS-RIS Channel
F=zeros(NL,NM);
theta_nmAOA= ((rand(1,Nf)*pi)-pi/2); 
theta_nmAOD=((rand(1,Nf)*pi)-pi/2); 
for iteaoa=1:1:length(theta_nmAOA)
    AoAF(:,iteaoa)=(exp(-1j*pi*sin(theta_nmAOA(iteaoa))*[0:NL-1])).';
    GenieaidedAoD(:,iteaoa)=(exp(-1j*pi*sin(theta_nmAOD(iteaoa))*[0:NM-1])).';
    F=F+sqrt(0.5)*(normrnd(0,1,1,1) + 1j*normrnd(0,1,1,1)).*(AoAF(:,iteaoa)* GenieaidedAoD(:,iteaoa)');%*exp(1j*rand(1,1) * 2 * pi); ()
end
F=sqrt(1/Nf).*F;

% RIS-User Channel
H=zeros(NL,K);
for itek=1:1:K
    theta_User=(rand(K,Nhk)*pi)-pi/2;%[10 30 ]';
    H(:,itek)=zeros(NL,1);
    for iteuaoa=1:1:Nhk
        AoDHk(:,iteuaoa,itek)=((exp(-1j*pi*sin(theta_User(itek,iteuaoa))*[0:NL-1])).');
        GenieaidedAoA(:,(iteuaoa-1)*length(theta_nmAOA)+1:(iteuaoa)*length(theta_nmAOA),itek)=diag(AoDHk(:,iteuaoa,itek)')*AoAF;
        H(:,itek)= H(:,itek)+sqrt(0.5)*(normrnd(0,1,1,1) + 1j*normrnd(0,1,1,1)).*AoDHk(:,iteuaoa,itek);
    end
    H(:,itek)= sqrt(1/(Nhk)).*H(:,itek);
end
end
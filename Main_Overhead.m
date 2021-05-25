clc
clear
close all
%% Code of the paper : £¨Figure 4£©
%J. Chen, Y.-C. Liang, H. V. Cheng, and W. Yu, ¡°Channel estimation for%reconfigurable intelligent surface aided multi-user mimo systems,¡± IEEEb Trans. Wireless Commun., 2021(accept).
%%
Channel_realization=50;% The numer of Monte Carlo Trials
ParameterB=[ 8 24  40   72   128 ];% The numer of subframes
System.NM=128; % The number of antennas at the BS
System.NL=128;  % The numer of reflective elements at the RIS
System.K=4;  % The numer of users
System.SNR=10.^(10/10); % Transmit power 10db
System.Gr=512; % AoA Dictionary Angular Resolutions
System.Gt=128; % AoD Dictionary Angular Resolutions
System.Nf=8; % The number of path between BS-RIS channel
System.Nhk=1;% The number of path between RIS-U_k channel
System.N_npilot=System.K; % The bumber of pilots of each user (parameter $T$)
System.Dr= ULAdictionaryGeneralize(System.NL,System.Gr);
System.Dt= ULAdictionaryGeneralize(System.NM,System.Gt);
Userpilot=dftmtx(System.K);  % The pilot sequence of all users
NM=System.NM; % The number of antennas at the BS
NL=System.NL;  % The numer of reflective elements at the RIS
K=System.K;  % The numer of users
SNR=System.SNR;
Gr=System.Gr;
Gt=System.Gt;
Nf=System.Nf;
Nhk=System.Nhk;
N_npilot=System.N_npilot;
Dr=System.Dr;
Dt=System.Dt;

NMSE_SMMVAll=zeros(Channel_realization,length(ParameterB));
NMSE_SMJCEAll=zeros(Channel_realization,length(ParameterB));
NMSE_LS_gen_subspace=zeros(Channel_realization,length(ParameterB));
for ite_channel=1:1:Channel_realization
    nmse_SMMV=zeros(1,length(ParameterB));
    nmse_SMJCE=zeros(1,length(ParameterB));
    nmse_LS_gen_subspace=zeros(1,length(ParameterB));
    [F,H,GenieaidedAoA,GenieaidedAoD]=ChannelGenralize(System) ;
    % BS-RIS channel, RIS user channel, exact AoAs and exact AoDs for S-genieaided LS
    
    for ite_index=1:1:length(ParameterB)
        N_bpilot=ParameterB(ite_index);
        % System.N_bpilot=N_bpilot;
        % V=RIS_SequenceOptimization(System,N_bpilot); %RIS Reflection Coefficient Optimization  
         LISpilot=sqrt(0.5)*(normrnd(0,1,NL,N_bpilot) + 1j*normrnd(0,1,NL,N_bpilot));% 
        while rank(LISpilot)<min(N_bpilot,NL)
            LISpilot=sqrt(0.5)*(normrnd(0,1,NL,N_bpilot) + 1j*normrnd(0,1,NL,N_bpilot)); 
        end
        LISpilot=LISpilot./abs(LISpilot);
        V=LISpilot;
        % Ignoring the direct linkr channelestimation
        % Received signal Model============================================
        Yk=zeros(NM,N_bpilot,K);
        noise=zeros(NM,N_bpilot,K);
        YkHermite=zeros(N_bpilot,NM,K);
        U=sqrt(0.5)*(normrnd(0,1,NM,N_npilot*N_bpilot) + 1j*normrnd(0,1,NM,N_npilot*N_bpilot));
        G=zeros(NL,NM,K);
        for itek=1:1:K
            G(:,:,itek)= diag(H(:,itek)')*F;
            for iteb=1:1:N_bpilot
                noise(:,iteb,itek)=(U(:,N_npilot*(iteb-1)+1:N_npilot*iteb)*Userpilot(:,itek))./(sqrt(SNR).*(N_npilot)) ;%.*sqrt(Nu*Nb).
                Yk(:,iteb,itek)= G(:,:,itek)'*V(:,iteb) +  noise(:,iteb,itek);
            end
            YkHermite(:,:,itek)=Yk(:,:,itek)'; %eq(16)
        end
        %% Proposed Algorithm =====================================
        [NMSE_SMJCE,NMSE_SMMV]=SMJCE( Yk, YkHermite,V,G,System,N_bpilot);
        nmse_SMJCE(1,ite_index)=NMSE_SMJCE; %The NMSE performance of the proposed algorithm.
        nmse_SMMV(1,ite_index)=NMSE_SMMV;
        
        %% Estimtated Error by Subspace Gen aided LS================
        [AoDsvdall,eigsvd,~]=svd(GenieaidedAoD);
        [~,length_indTRUE]=size(GenieaidedAoD);
        AoDsvd=AoDsvdall(:,1:length_indTRUE);
        z_bktildeinv_true=zeros(N_bpilot,length_indTRUE,K);
        G_hatall_gen= zeros(NL,NM,K);
        for itek=1:1:K
            [AoAsvdkall,~,~]=svd(GenieaidedAoA(:,:,itek));
            AoAsvdk= AoAsvdkall(:,1:Nhk*Nf);
            z_bktildeinv_true(:,:,itek)=YkHermite(:,:,itek)*AoDsvd*(AoDsvd'*AoDsvd)^(-1);
            G_hat=AoAsvdk*((V'*AoAsvdk)'*(V'*AoAsvdk))^(-1)*(V'*AoAsvdk)'*z_bktildeinv_true(:,:,itek)*AoDsvd';
            G_hatall_gen(:,:,itek)=G_hat;
            nmse_LS_gen_subspace(1,ite_index)=nmse_LS_gen_subspace(1,ite_index)+trace((G_hat-G(:,:,itek))'*(G_hat-G(:,:,itek)))/trace((G(:,:,itek))'*(G(:,:,itek)))/K;
        end 
    end
    
    NMSE_LS_gen_subspace(ite_channel,:)=nmse_LS_gen_subspace;
    NMSE_SMMVAll(ite_channel,:)=nmse_SMMV;
    NMSE_SMJCEAll(ite_channel,:)=nmse_SMJCE;
    fprintf('%dth channel is done...................\n',ite_channel)
end

NMSE_LS_gen_subspace_av=sum(NMSE_LS_gen_subspace(1:Channel_realization,:),1)./Channel_realization;
NMSE_SMJCE_av=sum(NMSE_SMJCEAll(1:Channel_realization,:),1)./Channel_realization;
NMSE_SMMV_av=sum(NMSE_SMMVAll(1:Channel_realization,:),1)./Channel_realization;



figure(2)
semilogy(ParameterB,NMSE_SMMV_av,'-bp','LineWidth',1.4,'MarkerSize',8);
hold on
semilogy(ParameterB,NMSE_SMJCE_av,'-kp','LineWidth',1.4,'MarkerSize',8);
semilogy(ParameterB,NMSE_LS_gen_subspace_av,'-rp','LineWidth',1.4,'MarkerSize',8);
legend('S-MMV','S-MJCE','S-Genie aided LS')
ylabel('NMSE')
xlabel('Training overhead {\itB}')
grid on
set(gca, 'FontName', 'Times New Roman')
box on
set(gca, 'FontName', 'Times New Roman')
set(gca,'FontSize',11); % ????????????????????????????????????????????????
set(get(gca,'XLabel'),'FontSize',10);%??????????8 point????5??
set(get(gca,'YLabel'),'FontSize',10);
grid on

%set(gca,'xtick',[8 24 40 56 72 88 104 128]);
%set(gca,'xticklabel',{'8','24','40','56','72','88','104','128'});
%xlim([8 128])
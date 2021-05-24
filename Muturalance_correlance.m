function  [Number_random,MeanC,varC]=Muturalance_correlance(Gr,V,Dr,xbins1)
ind=1:1:Gr;
        ind=(ind-1)*Gr+ind;
        Dic=V'*Dr; 
        for ite=1:1:Gr
        Dic(:,ite)=Dic(:,ite)./norm(Dic(:,ite));
        end 
        Grammatrix=(Dic'*Dic);
        vecGrammatrix=abs(vec(Grammatrix));
        vecGrammatrix(ind)=[];
      % yy2=dsort(vecGrammatrix);
        [Number_random,~]=hist(vecGrammatrix,xbins1);
        MeanC=mean(vecGrammatrix);
        varC=var(vecGrammatrix);
end
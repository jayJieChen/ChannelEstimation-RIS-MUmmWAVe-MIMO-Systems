function Dictionary=ULAdictionaryGeneralize(Dimension,NumResolu)
Dictionary=zeros(Dimension,NumResolu);  
for ite=1:1:NumResolu
    Dictionary(:,ite)=(sqrt(1/Dimension)*exp(-1j*2*pi*(-0.5+(ite-1)/NumResolu)*[0:Dimension-1])).';
end 
end
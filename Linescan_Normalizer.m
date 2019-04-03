% Normalize Linescan to Integrated intensity of 1st spectrum. 
Linescan = double(SI_Linescan2);
sizeL=size(Linescan);
Linescan_N=Linescan;
IntegratedZLP=zeros(sizeL(1),1);
for k = 1:sizeL(1)
    Linescan_N(k,:)=Linescan(k,:)*sum(Linescan(1,:))/sum(Linescan(k,:));
    IntegratedZLP(k,1)= sum(Linescan(k,:));
end
function threshmat=get_fdr_bullmore(mat)
%% Matrice is partial
z=.5.*log((1+mat)./(1-mat));

[~,p]=corrcoef(z);

p=triu(p);


p(find(p==0))=1;

threshmat=mat;
p=reshape(p,size(mat,1)*size(mat,1),1);
[rr,cc]=sort(p);

N=size(mat,1)*(size(mat,1)-1)/2; somme=0;
for i=1:N
    somme=somme+1/i;
end
index=[];
for i=1:length(rr)
    if(rr(i)>0.05*i/(N*somme))
      threshmat( cc(i))=0;
    end
end

 threshmat=threshmat+threshmat';
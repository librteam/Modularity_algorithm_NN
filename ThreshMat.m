function thresholded=ThreshMat(Matrix,percentage)

%% Matrix Size: nbchannels x nbchannels

numedges=nnz(Matrix(:,:,1));
neededEdges=round((percentage*numedges)/100);

if mod(neededEdges,2)~=0
    neededEdges=neededEdges+1;
end

% for i=1:size(Matrix,3)
clear b z; 
z=Matrix;
l=find(~z);
% z(l)=-inf;
[z,q]=sort(z(:));
notedge=size(Matrix,1)*size(Matrix,1)-numedges;
z(1:notedge+numedges-neededEdges)=0;
b(q)=z;
b(l)=0;
aa=reshape(b,size(Matrix,1),size(Matrix,1));
thresholded=aa;
% end
function [matrix, finalMatrix, start,fin]=code_seg(m,err)

[start1,fin1]=get_segmented(m,err);
matrix1=zeros(size(m));
for i=1:length(start1)
    matrix1(start1(i):fin1(i),start1(i):fin1(i))=1;
end

m2=m(size(m,1):-1:1,size(m,2):-1:1);

[s,f]=get_segmented(m2,err);
fin2=size(m,1)+1-s(length(s):-1:1);
start2=size(m,1)+1-f(length(f):-1:1);
matrix2=zeros(size(m));
for i=1:length(start2)
    matrix2(start2(i):fin2(i),start2(i):fin2(i))=1;
end
matrix = matrix1+matrix2;
matrix(matrix>0.5)=1;

[start,fin]=get_segmented(matrix,0.5);

% start=0;
% fin=0;

% dg=diag(matrix);

% start=[start1;start2];
% start=min(start);
% 
% fin=[fin1;fin2];
% fin=max(fin);


finalMatrix=zeros(size(m));
for i=1:length(start)
    finalMatrix(start(i):fin(i),start(i):fin(i))=1;
end

end

function [start,fin]=get_segmented(m,err)
% err=0.7;


start=[];
fin=[];
c=0;
for i=1:size(m,1)
    if c>1
        c=c-1;
        continue;
    end
    c=0;
    for j=i:size(m,1)
        if (sum(m(j,i:j))+sum(m(i:j,j))-m(j,j))/((j-i)*2+1)<err
            break;
        else
            c=c+1;
        end
    end
    if c>0
        start=[start, i];
        fin=[fin, i+c-1];
        %         i=i+c;
    end
    
end
% matrix=zeros(size(m));
% for i=1:length(start)
%     matrix(start(i):fin(i),start(i):fin(i))=1;
% end
end
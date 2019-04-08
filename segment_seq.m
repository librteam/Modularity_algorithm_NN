function [start,fin,f2]=segment_seq(conn_net,gamma)


%% Step1: Normalize the similarity matrix
conn_net1=(conn_net-min(min(conn_net)))/(max(max(conn_net))-min(min(conn_net)));
con_thresh=get_fdr_bullmore(conn_net1);

%% Step2: Median filter
i=0;
oldimage=con_thresh;
while(true)
    i=i+1
    newimage=medfilt2(oldimage);
    if(sum((sum(newimage-oldimage))~=0))
        oldimage=newimage;
    else
        break;
    end
end

%% Step3: Normalize the similarity matrix
newimage=(newimage-min(min(newimage)))/(max(max(newimage))-min(min(newimage)));


[f1,f2,start,fin]=code_seg(newimage,gamma);


function [conn_net,Modules]=Mod_AK_seq(plv_dyn,thresholdtype,thresholdperc)

%% Step1:Extract Dynamic modules
Ntw=size(plv_dyn,1);
Modules=[];
MatModules=[];
for tw=1:size(plv_dyn,1) %% over time windows
        threshplv_tw=squeeze(plv_dyn(tw,:,:));    
    if(strcmp(thresholdtype,'FDR'))
          % Apply statistical FDR threshold
            threshplv_tw=get_fdr_bullmore(threshplv_tw);
    else
        
        if (strcmp(thresholdtype,'prop'))
              % Apply proportional threshold
            threshplv_tw=ThreshMat(threshplv_tw,thresholdperc);
        end
    end
    
    trials=[];
    for iter=1:200
        M_Louvain = community_louvain(threshplv_tw);
        trials=[trials,M_Louvain];
    end
    [newM ,~, X_new3,~,~] = consensus_iterative(trials');
    Modules=[Modules;newM(1,:)];
    MatModules(tw,:,:)=X_new3;
end

%% Step2: Compute the similarity matrix between time windows
conn_net=zeros(Ntw,Ntw);
for i=1:Ntw
    i
    for j=i+1:Ntw
        [~ ,~,conn_net(i,j),~]=zrand(Modules(i,:),Modules(j,:));
    end
end
conn_net=conn_net+conn_net';

gamma=mean(mean(conn_net));

[start,fin,f2]=segment_seq(conn_net,gamma);


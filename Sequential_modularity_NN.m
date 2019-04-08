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

%% Step3: apply the segmentation algorithm 
[start,fin,f2]=segment_seq(conn_net,gamma);

%% Step4: Reject random modules and compute the final modules representing the dynamic networks

finalmodules=[];
final_allegiance=[];
for i=1:length(start)
    inst=[start(i):fin(i)];
    Mi=Modules(inst,:);
    [Module_i ,~, Allegiance_i,~,zrand] = consensus_iterative(Mi);
    finalmodules(i,:)=Module_i(1,:);
    final_allegiance(i,:,:)=Allegiance_i;
    %%rand is the distance from the null model computed from random permutations : 0==random; 1= no randomness
    rand(i)=zrand;
end

modindex=find(rand>0);
Imp_modules=finalmodules(modindex,:);
Imp_association=final_allegiance(modindex,:,:);

%% change the dynamic Module affiliation accordingly
Module_dynamic=zeros(1,201)
for i=1:length(start)
Module_dynamic(start(i):fin(i))=i;
end


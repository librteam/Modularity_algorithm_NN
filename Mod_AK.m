function [conn_net,Imp_modules, Imp_association, Module_dynamic]=Mod_AK(plv_dyn,thresholdtype,thresholdperc)

%% ---------Input---------- %%
% 1- plv_dyn: Dynamic plv matrix of size Ntw * Nnodes * Nnodes where Ntw is the
%number of time windows (FULL ADJ MATRIX)
% 2- threshold type: 'FDR' or 'prop' 
% 3- threshold percentage: a value ranging from 1 to 100% needed to apply
% the proportional threshold


%% --------Outputs----------%%
% 1-conn_net: Similarity matrix
% 2-Imp_modules: The modules representing the dynamic networks
% 3-Imp_association: The association matrices of Imp_modules
% 4-Module_dynamic: The module affiliation at each time window(from 0-->the number of
% Imp_modules)

%% Developed by Aya Kabbara 
%% ----------------------------------------------

Ntw=size(plv_dyn,1);

%% step1: Extract Dynamic modules
Modules=[];
MatModules=[];
%% time windows * ROIs * ROIs
for tw=1:size(plv_dyn,1) %% over time windows
    threshplv_tw=squeeze(plv_dyn(tw,:,:));
    
%     if(~issymmetric(plv_tw))
%     threshplv_tw=plv_tw+plv_tw';
%     end
    threshplv_tw(threshplv_tw<0)=0;
    
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
%imagesc(conn_net)

%% Step3: Get the partition of time windows into modules based on the similarity matrix
conn_net(conn_net<0)=0;
trials=[];
for iter=1:200
    M_Louvain = community_louvain(conn_net,1.1);
    trials=[trials,M_Louvain];
end
[Module_aff ,~, Allegiance_mat,~,~] = consensus_iterative(trials');
Module_aff=Module_aff(1,:);

%% Compute the association matrix of each temporal module
finalmodules=[];
final_allegiance=[];
for i=1:max(Module_aff)
    inst=find(Module_aff==i);
    Mi=Modules(inst,:);
    [Module_i ,~, Allegiance_i,~,rand1] = consensus_iterative(Mi);
    finalmodules(i,:)=Module_i(1,:);
    final_allegiance(i,:,:)=Allegiance_i;    
    %%rand is the distance from the null model computed from random permutations : 0==random; 1= no randomness
    rand(i)=rand1;
end

%% Step4: Reject random modules and compute the final modules representing the dynamic networks
modindex=find(rand>0);
Imp_modules=finalmodules(modindex,:);
Imp_association=final_allegiance(modindex,:,:);

%Change the dynamic Module affiliation accordingly
Module_dynamic=[];
for i=1:size(Module_aff,2)
    if(rand(Module_aff(i))==0)
        Module_dynamic(i)=0;
    else
        Module_dynamic(i)=find(modindex==Module_aff(i));
    end
end

% figure;
% imagesc(Module_dynamic)

% %% proportional occurrence
% for i=1:(max(Module_dynamic))
%     prop_occ(i)=max(max(squeeze(Imp_association(i,:,:))))/Ntw;
% end
% prop_occ
% rand

% %% Similarity between modules
% sim_modules=zeros(size(Imp_modules,1),size(Imp_modules,1));
% for i=1:size(Imp_modules,1)
%     for j=i+1:size(Imp_modules,1)
%         
%     [~,~,sim_modules(i,j),~]=zrand(Imp_modules(i,:),Imp_modules(j,:));
%     end
% end
% sim_modules=sim_modules+sim_modules';
% nx=size(sim_modules,1);ny=nx;
% figure;
% imagesc(sim_modules)
% colorbar
% title('Similarity matrix')
% set(gca,'xtick', linspace(0.5,nx+0.5,nx+1), 'ytick', linspace(0.5,ny+.5,ny+1));
% set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k');
% 
% %% transition matrix
% transi=Module_dynamic([1 diff(Module_dynamic)]~=0);
% taille=length(transi);
% Prob=zeros(max(Module_dynamic),max(Module_dynamic));
% 
% for i=1:max(Module_dynamic) %% t
%     for j=1:max(Module_dynamic) %% t+1
% %         if (i==j)
% %             Prob(i,j)=1;
% %         else
%             index_i=find(transi==i);
%             if(length(index_i)==0)
%                 Prob(i,j)=0;
%             else
%             index_i_post=index_i+1;
%             if(index_i_post(end)==taille+1)
%                 index_i_post=index_i_post(1:end-1);
%             end
%             
%             post=transi(index_i_post);
%             Prob(i,j)=length(find(post==j))/length(find(transi==i));
% %         end
%             end
%     end
% end
% nx=max(Module_dynamic);ny=max(Module_dynamic);
% image;
% imagesc(Prob)
% colorbar
% title('Transition Matrix')
% set(gca,'xtick', linspace(0.5,nx+0.5,nx+1), 'ytick', linspace(0.5,ny+.5,ny+1));
% set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k');
% 

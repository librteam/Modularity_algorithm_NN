%% load simulated dynamic matrices 

load(S_noise);

%% Generate categorical MS
[similarity_mat,Imp_modules, Imp_association, dynamic_affiliation]=categorical_modularity_NN(plv_dyn,'prop',90)

% Display the dynamic affiliations
imagesc(dynamic_affiliation);

%% Generate consecutive MS
[conn_net,Imp_modules, Imp_association, Module_dynamic]=sequential_modularity_NN((plv_dyn,'prop',90)

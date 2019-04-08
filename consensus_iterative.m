function [S2, Q2, X_new3, qpc, disrand] = consensus_iterative(C)
%CONSENSUS_ITERATIVE     Construct a consensus (representative) partition
%using the iterative thresholding procedure
%
%   [S2 Q2 X_new3 qpc] = CONSENSUS_ITERATIVE(C) identifies a single
%   representative partition from a set of C partitions, based on
%   statistical testing in comparison to a null model. A thresholded nodal
%   association matrix is obtained by subtracting a random nodal
%   association matrix (null model) from the original matrix. The
%   representative partition is then obtained by using a Generalized
%   Louvain algorithm with the thresholded nodal association matrix.
%
%   NOTE: This code requires genlouvain.m to be on the MATLAB path
%
%   Inputs:     C,      pxn matrix of community assignments where p is the
%                       number of optimizations and n the number of nodes
%
%   Outputs:    S2,     pxn matrix of new community assignments
%               Q2,     associated modularity value
%               X_new3, thresholded nodal association matrix
%               qpc,    quality of the consensus (lower == better)
%               
%
%   Bassett, D. S., Porter, M. A., Wymbs, N. F., Grafton, S. T., Carlson,
%   J. M., & Mucha, P. J. (2013). Robust detection of dynamic community
%   structure in networks. Chaos: An Interdisciplinary Journal of Nonlinear
%   Science, 23(1), 013142.

npart = numel(C(:,1)); % number of partitions
m = numel(C(1,:)); % size of the network

% initialize
C_rand3 = zeros(size(C)); % permuted version of C
X = zeros(m,m); % Nodal association matrix for C
X_rand3 = X; % Random nodal association matrix for C_rand3

%% NODAL ASSOCIATION MATRIX

% try a random permutation approach
for i = 1:npart;
    pr = randperm(m);
    C_rand3(i,:) = C(i,pr); % C_rand3 is the same as C, but with each row permuted
end

% Calculate the nodal association matrices X and
% X_rand3
for i = 1:npart;
    ii=i
    tic
    for k = 1:m
        for p = 1:m;
            % element (i,j) indicate the number of times node i and node j
            % have been assigned to the same community
            if isequal(C(i,k),C(i,p))
                X(k,p) = X(k,p) + 1;
            else
                X(k,p) = X(k,p) + 0;
            end
            
            % element (i,j) indicate the number of times node i and node j
            % are expected to be assigned to the same community by chance
            if isequal(C_rand3(i,k),C_rand3(i,p))
                X_rand3(k,p) = X_rand3(k,p)+1;
            else
                X_rand3(k,p) = X_rand3(k,p)+ 0;
            end
        end
    end
    toc
end

%% THRESHOLDING
% % keep only associated assignments that occur more often than expected in
% % the random data

before=length(find(tril(X,-1)));
  X_new3 = zeros(m,m);
% X_new3=X;
   X_new3(X>max(max(triu(X_rand3,1)))) = X(X>max(max(triu(X_rand3,1))));
   %X_new3(X>mean(mean(triu(X_rand3,1)))) = X(X>mean(mean(triu(X_rand3,1))));
after=length(find(tril(X_new3,-1)));

disrand=after/before;
%% GENERATE THE REPRESENTATIVE PARTITION
% recompute optimal partition on this new matrix of kept community
% association assignments
for i = 1:npart;
    ii=i
   [S2(i,:) Q2(i)] = multislice_static_unsigned(X_new3,1);

end

% define the quality of the consensus
qpc = sum(sum(abs(diff(S2))));

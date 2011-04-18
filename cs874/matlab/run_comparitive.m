function tids = run_comparative(X, Y, ile, ss)
%For the standard run, input both x: the gene expression values per tissue
%with genes indexed in the first axis and tissues on the second. As well as
%Y which has tissues on the first axis.

[U,S,V] = svd(X);
XV = X*V; XV=XV';

pca_sims=0;
mve_sims=0;
if pca_sims == 1
    sims = similarity(XV);
elseif mve_sims == 1
    sims = similarity(Y);
else
    sims = similarity(X);
end

if nargin < 3
  ile = 20;
end

if nargin < 4
    ss  = prctile(sims(:),ile);
end
   
size(sims)
tids  = ap(sims,ss);
show_results(XV,Y,tids);







function show_results(X,Y,clusters)
  %Show the results of affinity propagation in two dimensions in both the MVE
  %projection and the PCA projection computed with the SVD.
    cols = rand(max(clusters),3);
    subplot(2,1,1);
    scatter(Y(1,:),Y(2,:),100,cols(clusters,:), 'filled');
    
    subplot(2,1,2);
    scatter(X(1,:),X(2,:),100,cols(clusters,:), 'filled');

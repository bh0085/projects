
function [L] = cluster_match(ids1,ids2)

n = length(ids1); % total number of profiles

% unique cluster ids for clustering 1
N=max(ids1); Z=zeros(1,N); Z(ids1) = 1; C1 = find(Z); 

% unique cluster ids for clustering 2
N=max(ids2); Z=zeros(1,N); Z(ids2) = 1; C2 = find(Z); 

% match clusters pairwise 
L = []; 
for i=1:length(C1), 
   for j=1:length(C2),
      I1 = find(ids1==C1(i)); 
      I2 = find(ids2==C2(j)); 
      Z = zeros(1,n); Z(I1) = 1; 
      k = sum(Z(I2)); % overlap 
      n1= length(I1); % cluster 1 size
      n2= length(I2); % cluster 2 size
      L=[L;C1(i),C2(j),n1,n2,k]; 
   end;
end;
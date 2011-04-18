% function [S] = similarity(X)
% produces a similarity matrix containing negative squared distances between rows of X 

function [S] = similarity(X)

[n,d] = size(X); 
X2 = sum(X.^2,2); 
S = -repmat(X2,1,n)+2*X*X'-repmat(X2',n,1); 
S = S - diag(diag(S));


% function [bic] = bic_score(X,ids)
% evaluates the BIC score for a reconstructed mixture model corresponding
% to exemplar assignments in ids. Rows of X are interpreted as points.

function [bic] = bic_score(X,ids)

[n,d] = size(X); 

% overall variance
var_est = sum(sum( (X-X(ids,:)).^2 ))/(n*d); 

% mixing proportions
pi_est = zeros(n,1); 
for i=1:n, pi_est(ids(i)) = pi_est(ids(i)) + 1; end;
ex = find(pi_est); pi_est = pi_est(ex)/sum(pi_est); 
k = length(ex);

LL = 0; % log-likelihood 
for i=1:n,
    % log-likelihood of X(i,:) according to a mixture of Gaussians
    % model with means X(ex,:), overall variance var_est, and mixing
    % proportions pi_est 
    
    
    %gpdf = normpdf(repmat(X(i,:),k,1), X(ex,:),repmat(var_est,k,size(X,2)))
    %gauss1 = prod(gpdf,2);
    
    nrm = (2 * var_est*pi).^(-1*d/2);
    gauss2 = nrm .* exp(-.5 * sum((repmat(X(i,:), k,1) - X(ex,:)).^2,2)./ var_est);
    loglik=log(sum(pi_est .* prod(gauss2,2)));
    if isfinite(loglik) == 0; 
            bic = -inf;
            return
    end
    LL = LL + loglik; 
end; 

bic = LL - 0.5*(k*d+k-1+1)*log(n);


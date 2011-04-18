function [ids,ex] = ap(S,self_sim,constraint_list)

n = size(S,1); 
if (nargin<2), 
   for i=1:n, S(i,i) = median(S(:,i)); end; 
else
   for i=1:n, S(i,i) = self_sim; end;  
end;

% connected components (points assigned independently)
cc = cell(n,1); for i=1:n, cc{i} = i; end;
if (nargin==3), % generate connected components from constraint_list
   z = zeros(n,1); for i=1:n, z(i)=i; end;
   for i=1:size(constraint_list,1),
       zt = z(constraint_list(i,:)); zmin = min(zt); 
       z(find(z==zt(1)))=zmin; 
       z(find(z==zt(2)))=zmin;
   end;
   for i=1:n, cc{i} = find(z==i); end;
end;

lm = 0.25; % damping factor 

R = zeros(n,n); % responsibilities
A = zeros(n,n); % availabilities 
C = zeros(n,n); % biases from additional constraints

cont = 1; iter = 1; max_probs = diag(R+S+C+A); 
while (cont), 
    old_max_probs = max_probs; 

    % udpate responsibilities 
    AS = A+C+S; [Y,I]=max(AS,[],2);
    for i=1:n AS(i,I(i))=-Inf; end;
    [Y2,I2]=max(AS,[],2);
    Y=repmat(Y,[1,n]); for i=1:n Y(i,I(i))=Y2(i); end;
    R=(1-lm)*R-lm*Y; 
    
    % update availabilities 
    Rp=max(R+C+S,0); for i=1:n, Rp(i,i)=R(i,i)+S(i,i)+C(i,i); end;
    Anew=repmat(sum(Rp,1),[n,1])-Rp;
    dA=diag(Anew); Anew=min(Anew,0); for i=1:n, Anew(i,i)=dA(i); end;
    A=(1-lm)*A+lm*Anew;

    % enforce the same assignment for each point in a connected component
    F = R+S+A; 
    for i=1:n, 
       if (length(cc{i})>1), 
          Cnew = sum(F(cc{i},:),1); % pool the info from points within cc
          for j=1:length(cc{i}),
             ii=cc{i}(j); 
             C(ii,:) = (1-lm)*C(ii,:)+lm*(Cnew-F(ii,:));
          end;
       end;
    end;
    
    % stopping criterion
    max_probs = diag(R+S+C+A);
    cont = max(abs(max_probs-old_max_probs))>1e-6; 
    if (iter>10*n), cont = 0; disp('max number of iterations reached'); end; 
    iter = iter+1;
end;
 
[tmp,ex] = sort(max_probs,'descend');
k = max(sum(tmp>0),1);
ex = ex(1:k); % k exemplars

ids = zeros(n,1); % assignments to closest exemplars
for i=1:n,
   % find the best common assignment within the connected component
   [tmp,ccid] = max(sum(S(cc{i},ex),1),[],2); 
   % represent as individual assignments
   ids(cc{i}) = ccid; 
end;
ids(ex) = 1:length(ex);
ids = ex(ids); 

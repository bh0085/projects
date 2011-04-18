% function [] = plotpoints(X,ids,dim)
% plots two or three dimensional projections of points represented by rows of X
% connects each point to its exemplar id if provided

function [] = plotpoints(X,ids,dim)
    
[n,d] = size(X); 
if (nargin<3), dim = 2; end;
dim = min([dim,3,d]); 

if (d>2), 
   [U,S,V] = svd(X);
   X = X*V(:,1:dim);
end;

if (dim==3), 
   plot3(X(:,1),X(:,2),X(:,3),'o'); 
else
   plot(X(:,1),X(:,2),'o'); 
end;
hold on;

if (nargin>1), 
   for i=1:n,
      if (dim==3), 
         plot3([X(i,1),X(ids(i),1)],[X(i,2),X(ids(i),2)],[X(i,3),X(ids(i),3)],'r');
      else
         plot([X(i,1),X(ids(i),1)],[X(i,2),X(ids(i),2)],'r');
      end;
   end;
end;

hold off;


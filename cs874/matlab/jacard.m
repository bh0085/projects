% function [J] = jacard(ids1,ids2)
% evaluates the jacard index between two clusterings represented
% by vectors ids1 and ids2 where ids1(i) = cluster id of point i.

function [J] = jacard(ids1,ids2)

n = length(ids1); 

N11 = 0; N10 = 0; N01 = 0; 
for i=1:n,
   for j=i+1:n,
      if (ids1(i)==ids1(j) && ids2(i)==ids2(j)), N11 = N11 + 1; end;
      if (ids1(i)==ids1(j) && ids2(i)~=ids2(j)), N10 = N10 + 1; end;
      if (ids1(i)~=ids1(j) && ids2(i)==ids2(j)), N01 = N01 + 1; end;
   end;
end;
J = N11/(N11+N10+N01); 
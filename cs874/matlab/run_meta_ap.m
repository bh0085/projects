function [ks,bs, ss] = run_meta_ap(X)

sims = similarity(X);
ss = -1*logspace(log10(-prctile(sims(:),90)),log10(-min(sims(:))),30); 
ss = ss';

bs=zeros(size(ss,1),1);
ks=zeros(size(ss,1),1);


for i = 1:size(ss,1); 
    ids = ap(sims,ss(i));
    bs(i) = bic_score(X,ids);
    exemplars = zeros(size(ids,1),1);
    for j = 1:size(ids,1);
        exemplars(ids(j)) = exemplars(ids(j))+ 1;
    end
    ex_ids = find(exemplars);
    ks(i) = size(ex_ids,1);
end


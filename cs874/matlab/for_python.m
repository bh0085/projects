function for_python(X,Y,Y8,Y12)

[U,S,V] = svd(X);
XV = X*V; XV=XV';
bestss = -757.7858;

dorun_meta = 0;
if dorun_meta
    [ks,bs,ss] = run_meta_ap(X);
    [srt, idxs]= sort(bs, 'descend');
    ssval = ss(idxs(0));
else
    ssval = bestss;
    best_ile= 45.;
end

outp = {};

X = X;
XV= XV';
Y = Y';
Y8= Y8';
Y12= Y12';

data = open('miRNA.mat');

sims = similarity(X); ss = prctile(sims(:),best_ile);
ids_X = ap(sims,ss);

iter_doms = 0;
if iter_doms ==1 
    sims = similarity(XV(:,1:2)); ss = prctile(sims(:),best_ile);
    ids_XV = ap(sims,ss);

    sims = similarity(Y(:,1:2)); ss = prctile(sims(:),best_ile);
    ids_Y = ap(sims,ss);

    sims = similarity(Y8(:,1:2)); ss = prctile(sims(:),best_ile);
    ids_Y8 = ap(sims,ss);

    sims = similarity(Y12(:,1:2)); ss = prctile(sims(:),best_ile);
    ids_Y12 = ap(sims,ss);
else
    constraints = data.constraint_list;
    %lens = floor(linspace(5, size(constraints,1), 20)); lens = lens'
    lens = [5,6,7,10,15, 25, 50, 100, 150, 200, 300, 500, 700]'
    all_ids = zeros(size(ids_X,1),size(lens,1));
    for i = 1: size(lens,1)
        all_ids(:,i) = ap(sims,ssval, constraints(1:lens(i),:));
    end
end
   

expression = data.expression;
tissue_description = data.tissue_description;
tissue_category = data.tissue_category;
constraint_list = data.constraint_list;

save('prob2.mat')


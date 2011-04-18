function L = p22(ids)
  best_k = 19;
  best_ss = -757.7858;
  data = open('miRNA.mat');
  tcats = data.tissue_category;
  
  n = size(ids,1)
  ex_ids = zeros(n);
  for i=1:n; ex_ids(ids(i)) = 1; end
  ex_list = find(ex_ids);  
  ids_inlist = zeros(n,1)
  for i=1:n; ids_inlist(i)=find(ex_list==ids(i)); end
  L = cluster_match(tcats,ids_inlist);
 
  
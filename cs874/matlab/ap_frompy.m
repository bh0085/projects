function ap_frompy(mat_inp, mat_out)
data = open(mat_inp)

sims = data.similarities
ss   = data.self_similarity
inds = ap(sims, ss)

out_struct = struct()
out_struct.inds = inds
out_struct.ss = ss
out_struct.metric = data.metric

save(mat_out, 'out_struct')
exit
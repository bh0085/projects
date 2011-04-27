function ap_frompy(mat_inp, mat_out)
data = open(mat_inp)

out_struct = struct()
out_struct.ss = data.self_similarity
out_struct.metric = data.metric
out_struct.inds = []
save(mat_out, 'out_struct')

sims = data.similarities
ss   = data.self_similarity
inds = ap(sims, ss)
out_struct.inds = ap(sims, ss)
out_struct = struct()
out_struct.ss = ss
out_struct.metric = data.metric
save(mat_out, 'out_struct')
exit
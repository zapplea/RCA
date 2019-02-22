in_fname_ls = ['/home/yibing/Documents/data/bioinfo/RCA/RIPC_5TF_target.txt',
               '/home/yibing/Documents/data/bioinfo/RCA/RIPC_5mir_target.txt']
out_fname_ls = ['/home/yibing/Documents/data/bioinfo/RCA/5TF_input_binding.txt',
                '/home/yibing/Documents/data/bioinfo/RCA/5mir_input_binding.txt']

for i in range(len(in_fname_ls)):
    in_fname = in_fname_ls[i]
    out_fname = out_fname_ls[i]
    dic = {}
    with open(in_fname,'r') as f:
        for line in f:
            gene,tf = line.split('\t')
            tf = tf.replace('\n','')
            if gene in dic:
                dic[gene].append(tf)
            else:
                dic[gene]=['na',tf]
    for gene in dic:
        print('%s : %d'%(gene,len(dic[gene])))
    with open(out_fname,'w') as f:
        for gene in dic:
            tf_ls = dic[gene]
            tf_ls.insert(0,gene)
            s = '\t'.join(tf_ls)
            f.write(s+'\n')
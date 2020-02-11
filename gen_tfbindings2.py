# in_fname_ls = ['/home/yibing/Documents/data/bioinfo/RCA_new/RIPC_5TF_target_new.txt',
#                '/home/yibing/Documents/data/bioinfo/RCA_new/RIPC_5mir_target_new.txt',
#                '/home/yibing/Documents/data/bioinfo/RCA_new/RIPC_5TF_lncR.txt']
# out_fname_ls = ['/home/yibing/Documents/data/bioinfo/RCA_new/5TF_input_binding.txt',
#                 '/home/yibing/Documents/data/bioinfo/RCA_new/5mir_input_binding.txt',
#                 '/home/yibing/Documents/data/bioinfo/RCA_new/5TF_lncR_input_binding.txt']

in_fname_ls = ['/home/yibing/Documents/data/bioinfo/RCA_2020_2_10/untitled']
out_fname_ls = ['/home/yibing/Documents/data/bioinfo/RCA_2020_2_10/lncR_TF_input_binding.txt']

for i in range(len(in_fname_ls)):
    in_fname = in_fname_ls[i]
    out_fname = out_fname_ls[i]
    dic = {}
    with open(in_fname,'r') as f:
        for line in f:
            gene = line.split('\t')[0]
            tf = line.split('\t')[1]
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
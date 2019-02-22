in_fname = 'input_TF_binding.txt'

with open(in_fname,'r') as f:
    for line in f:
        print(repr(line))
        print('\n')

in_fname = '/home/yibing/Documents/data/bioinfo/RCA/input_TF_binding.txt'

with open(in_fname,'r') as f:
    for line in f:
        ls = line.split()
        print('%s : %d'%(ls[0],len(ls[1:])))
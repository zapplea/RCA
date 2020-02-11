from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np

def read_data(fname):
    with open(fname) as f:
        count =0
        rows = []
        data = []
        for line in f:
            line = line.replace('\n','')
            ls = line.split('\t')
            if count ==0:
                count+=1
                columns = ls[2:]
                continue
            value = list(map(float,ls[2:]))
            data.append(value)
            # TODO: convert _ to '\t'
            rows.append(ls[0]+'\t'+ls[1])

        df = pd.DataFrame(data,index=rows,columns=columns)
        return df

def calculate_pval(true_df,random_df_ls):
    column_name = true_df.columns.values
    row_name = true_df.index.values
    true_values = true_df.values
    cmp_result_ls = []
    for rdf in random_df_ls:
        rdf = rdf[column_name].loc[row_name]
        random_values = rdf.values
        condition = np.greater(true_values,random_values)
        cmp_result = np.where(condition,np.ones_like(random_values).astype('float32'),np.zeros_like(random_values).astype('float32'))
        cmp_result_ls.append(cmp_result)
    temp = np.sum(cmp_result_ls,axis=0)
    p_pos = 1-np.true_divide(temp,1000)
    p_neg = np.true_divide(temp,1000)
    condition = np.greater_equal(true_values,np.zeros_like(true_values))
    p_val = np.where(condition,p_pos,p_neg)
    return p_val,column_name.tolist(),row_name.tolist()

def write(fname,p_val,columns,rows):
    column = ['UNIQID\tNAME',]
    column.extend(columns)
    columns=column
    with open(fname,'w+') as f:
        f.write('\t'.join(columns)+'\n')
        p_val = p_val.tolist()
        for i in range(len(rows)):
            out_ls = [rows[i],]
            out_ls.extend(p_val[i])
            out_ls = list(map(str,out_ls))
            f.write('\t'.join(out_ls)+'\n')

if __name__ == '__main__':
    rootpath = '/home/yibing/Documents/data/bioinfo/RCAresult'
    true_fnames = [f for f in listdir(rootpath) if isfile(join(rootpath, f))]
    for fname in true_fnames:
        true_df = read_data(join(rootpath,fname))
        random_fdname = 'RCA_'+fname.split('.')[0]
        random_rootpath = join(rootpath,random_fdname)
        random_fnames = []
        for f in listdir(random_rootpath):
            if isfile(join(random_rootpath, f)):
                if f.split('.')[-1] == 'dat':
                    random_fnames.append(f)
        random_df_ls = []
        for rfname in random_fnames:
            rdf = read_data(join(random_rootpath,rfname))
            random_df_ls.append(rdf)
        p_value,columns,rows = calculate_pval(true_df,random_df_ls)
        write(join(rootpath,'p_'+fname.split('.')[0]+'.dat'),p_value,columns,rows)
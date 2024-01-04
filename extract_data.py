import pandas as pd
DEF_SZ: int = 10000

df_iter = pd.read_csv('./raw_data/raw_data_2.csv', header = [0, 1], chunksize=DEF_SZ, iterator=True)

out = open('./processed_data/processed_data_2.csv', 'w')

binding = dict()


for iter_num, chunk in enumerate(df_iter, 1):
    print(iter_num)
    if iter_num > 418:
        break
    for i in range((iter_num - 1) * DEF_SZ, iter_num * DEF_SZ, 1):
        #print(chunk[('Epitope', 'Name')][i], file = out)
        #print(chunk[('MHC Restriction', 'Name')][i], file = out)
        classs = chunk[('MHC Restriction', 'Class')][i]
        peptide = chunk[('Epitope', 'Name')][i]
        allele = chunk[('MHC Restriction', 'Name')][i]
        #print(peptide)
        #print(allele)
        if classs == "II" or len(peptide) != 9:
            continue
        if allele in binding:
            binding[allele].append(peptide)
        else:
            binding[allele] = list()
            binding[allele].append(peptide)
        #print(binding[allele])


for key in binding:
    if len(binding[key]) < 50:
        continue
    print(key, file = out)
    for i in binding[key]:
        print(i, end = ' ', file = out)
    print('\n', file = out)
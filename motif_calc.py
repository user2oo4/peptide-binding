
import linecache

output_path_motif = 'motif_data/motif_data_2.txt'
output_path_dist = 'motif_data/dist_data_2.txt'
output_path_debug = 'output.txt'
input_path = 'processed_data/processed_data_2.csv'
out_motif = open(output_path_motif, 'w')
out_debug = open(output_path_debug, 'w')

binding = dict()
motif = dict()
dist = dict()

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'W', 'Y', 'V']

def count_lines(file_path):
    with open(file_path, 'r') as file:
        line_count = sum(1 for line in file)
    return line_count

linesCnt = count_lines(input_path)
for i in range (1, linesCnt, 3):
    allele = linecache.getline(input_path, i)
    raw_peptides = linecache.getline(input_path, i + 1)
    #print(allele, file = out)
    if allele.find('*') == -1 and len(allele) != 6:
        continue
    cur = ""
    peptides = list()
    for c in raw_peptides:
        if c == ' ':
            peptides.append(cur)
            cur = ""
        else:
            cur += c
    #print(peptides, file = out)
    binding[allele] = peptides

for allele in binding:
    print(allele, file = out_motif)
    peptides = binding[allele]
    cnt = list(dict())
    sum = len(peptides)
    for i in range(0, 9):
        cnt.append(dict())
    for i in range(0, 9):
        for aa in amino_acids:
            cnt[i][aa] = 0
    for peptide in peptides:
        if peptide.find('X') != -1:
            continue
        for i in range(0, 9):
            try:
                cnt[i][peptide[i]] += 1
            except:
                print(peptide)
    for i in range(0, 9):
        for aa in cnt[i]:
            cnt[i][aa] /= sum
            print(cnt[i][aa], file = out_motif)
    motif[allele] = cnt


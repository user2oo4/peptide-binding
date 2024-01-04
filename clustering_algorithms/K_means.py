import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances
import os
from scipy.stats import entropy

def jensen_shannon_divergence(p, q):
    """Calculate Jensen-Shannon Divergence."""
    m = 0.5 * (p + q)
    return 0.5 * (entropy(p, m) + entropy(q, m))

def jensen_shannon_variance(distributions):
    """Calculate Jensen-Shannon Variance."""
    num_distributions = len(distributions)
    variance_matrix = np.zeros((num_distributions, num_distributions))

    for i in range(num_distributions):
        for j in range(i, num_distributions):
            jsd = jensen_shannon_divergence(distributions[i], distributions[j])
            variance_matrix[i, j] = jsd
            variance_matrix[j, i] = jsd

    variance = np.var(variance_matrix)
    return variance

supertypes = list()
os.environ["OMP_NUM_THREADS"] = "1"

print('num_groups??')
num_groups = int(input())
print('distance_type?? 0: Euclid, 1: Pearson, 2: JSD')
distance_type = int(input())

input_path = '../motif_data/motif_data_2.txt'
output_path = '../clustering_data/K_means_' + str(num_groups) + '_' + str(distance_type) + '.txt'
out_file = open(output_path, 'w')
in_file = open(input_path, 'r')

ord = dict() # map alleles with number
inv = dict() # inverse ord


num_alleles = 176 # number of alleles
num_positions = 9 # number of positions
num_aa = 20 # number of amino acides

data_dict = dict()


for i in range (0, num_alleles):
    str = next(in_file).strip()
    while str == '':
        str = next(in_file).strip()
    allele = str
    print(allele)
    ord[i] = allele
    inv[allele] = i
    for j in range (0, num_positions):
        for k in range(0, num_aa):
            str = next(in_file).strip()
            while str == '':
                str = next(in_file).strip()
            data_dict[i, j, k] = float(str)
# Convert the dictionary to a list of lists
data_list = [[[data_dict[i, j, k] for k in range(num_aa)] for j in range(num_positions)] for i in range(num_alleles)]

#print(data_list)
#exit(0)

# Compute the pairwise distances
distances = np.zeros((num_alleles, num_alleles))
if distance_type == 0:
    for i in range(0, num_alleles):
        for j in range(i + 1, num_alleles):
            for pos in range(0, num_positions):
                res = 0
                for aa in range(0, num_aa):
                    res += (data_list[i][pos][aa] - data_list[j][pos][aa]) * (data_list[i][pos][aa] - data_list[j][pos][aa])
                distances[i][j] += res ** 0.5
            distances[j][i] = distances[i][j]
elif distance_type == 1:
    for i in range(0, num_alleles):
        for j in range(i + 1, num_alleles):
            arr1 = np.array(data_list[i])
            arr2 = np.array(data_list[j])
            distances[i][j] = 1 - np.corrcoef(arr1, arr2)[0, 1]
            distances[j][i] = distances[i][j]
else:
    for i in range(0, num_alleles):
        for j in range(i + 1, num_alleles):
            for pos in range(0, 9):
                arr1 = np.array(data_list[i][pos])
                arr2 = np.array(data_list[j][pos])
                distances[i][j] += jensen_shannon_variance([arr1, arr2])
            distances[j][i] = distances[i][j]

#exit(0)
# Use k-means with precomputed distances
kmeans = KMeans(n_clusters = num_groups, init='random', n_init = 1)
kmeans.fit_predict(distances)

# Get cluster assignments and centroids
cluster_assignments = kmeans.labels_
centroids = kmeans.cluster_centers_

print(len(cluster_assignments))
# Print the results
for i in range (0, num_groups):
    supertypes.append(list())

for i in range(0, len(cluster_assignments)):
    print(ord[i], cluster_assignments[i])
    supertypes[cluster_assignments[i]].append(i)

for i in range(0, num_groups):
    print('Supertype', i, sep = ' ', file = out_file)
    for j in supertypes[i]:
        allele = ord[j]
        print(allele, file = out_file)
    print(file = out_file)

min_dist = list()
avg_dist = list()

for g1 in range(0, num_groups):
    min_dist.append(list())
    avg_dist.append(list())
    for g2 in range(0, num_groups):
        avg_dist[g1].append(0)
        min_dist[g1].append(100)
    for g2 in range(g1, num_groups):
        for a1 in supertypes[g1]:
            for a2 in supertypes[g2]:
                min_dist[g1][g2] = min(min_dist[g1][g2], distances[a1][a2])
                avg_dist[g1][g2] += distances[a1][a2]
        avg_dist[g1][g2] /= (len(supertypes[g1]) * len(supertypes[g2]))

        #print('min distance', g1, g2, '=', min_dist[g1][g2], sep = ' ')
        print('avg distance', g1, g2, '=', avg_dist[g1][g2], sep = ' ', file = out_file)


#print("Cluster Assignments:", cluster_assignments)
#print("Centroids:", centroids)

# Distance = 1 - Pearson correlation
# Jensen-Shannon distribution
# Only use data from mass-spec experiments
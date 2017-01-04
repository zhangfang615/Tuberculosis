from __future__ import division
def get_cluster_truth(clusters_truth_file):
    cluster_truth = {}
    line = clusters_truth_file.readline().strip()
    i = 0
    while line:
        patients = line.split(" ")
        for patient in patients:
            cluster_truth[patient]=i
        i += 1
        line = clusters_truth_file.readline().strip()
    return cluster_truth

def get_cluster_count(cluster_set,cluster_truth):
    cluster_counts = []
    for cluster in cluster_set:
        cluster_count = {}
        for patient in cluster:
            cluster = cluster_truth[patient]
            if cluster in cluster_count.keys():
                cluster_count[cluster] += 1
            else:
                cluster_count[cluster] = 1
        cluster_counts.append(cluster_count)
    return cluster_counts


def combination(number):
    if number > 1:
        return number*(number-1)/2
    else:
        return 0

def claculate_TPTN(cluster_numbers):
    TP = 0
    TN = 0
    FP = 0
    total = 0
    for i in range(len(cluster_numbers)):
        total_local = 0
        for cluster in cluster_numbers[i].keys():
            total += cluster_numbers[i][cluster]
            total_local += cluster_numbers[i][cluster]
            TP += combination(cluster_numbers[i][cluster])
            tn = 0
            for j in range(len(cluster_numbers)):
                if i != j:
                    for c in cluster_numbers[j].keys():
                        if c != cluster:
                            tn += cluster_numbers[j][c]
            TN += cluster_numbers[i][cluster]*tn
        FP += combination(total_local)
    return TP,TN,FP,total

resi_distance_matrix_file = file("./output_new/distance_unresi.txt")
line = resi_distance_matrix_file.readline().strip()
resi_index = line.split("\t")[1:]
size = len(resi_index)
resi_distance_matrix = []
line = resi_distance_matrix_file.readline().strip()
while line:
    distances = line.split("\t")[1:]
    resi_distance_matrix.append(distances)
    line = resi_distance_matrix_file.readline().strip()

resi_distance_matrix_file.close()

cluster_set = []
global_patient_set = set()
for i in range(size):
    patient = resi_index[i]
    cluster = set()
    if patient not in global_patient_set:
        cluster.add(patient)
        global_patient_set.add(patient)
        for j in range(i + 1, size):
            distance = int(resi_distance_matrix[i][j])
            if distance <75 and distance > -1 and resi_index[j] not in global_patient_set:
                cluster.add(resi_index[j])
                global_patient_set.add(resi_index[j])
    if len(cluster) > 0:
        cluster_set.append(cluster)


clusters_truth_file = file("./output_new/0.26_100_0.01_0.018_0.05_0.75_unresi.txt")
cluster_truth = get_cluster_truth(clusters_truth_file)

cluster_numbers = get_cluster_count(cluster_set,cluster_truth)

# cluster_numbers = []
# cluster1 = {}
# cluster1[1] = 5
# cluster1[2] = 1
# cluster_numbers.append(cluster1)
#
# cluster1 = {}
# cluster1[1] = 1
# cluster1[2] = 4
# cluster1[3] = 1
# cluster_numbers.append(cluster1)
#
# cluster1 = {}
# cluster1[1] = 2
# cluster1[3] = 3
# cluster_numbers.append(cluster1)

TP,TN,FP_TP,total = claculate_TPTN(cluster_numbers)
TN /= 2
print TP
print TN
print FP_TP
print total

RI = (TP+TN)/combination(total)
p = TP/FP_TP
r = TP/(combination(total)-FP_TP-TN+TP)
F = 2*p*r/(p+r)
print RI
print p
print r
print F
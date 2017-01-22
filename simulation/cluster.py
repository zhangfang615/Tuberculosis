from __future__ import division
from random import shuffle
import random
import pandas

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

def check_convergence(resi_index,patient_cluster,patient,global_patient_set,resi_distance_matrix):
    removed = []
    for patient_i in patient_cluster:
        if patient_i != patient:
            i = resi_index.index(patient_i)
            minimun = min(resi_distance_matrix[i])
            smallest = resi_distance_matrix[i].index(minimun)
            s = resi_index[smallest]
            if s not in patient_cluster:
                removed.append(patient_i)

    for patient in removed:
        patient_cluster.remove(patient)
        global_patient_set.remove(patient)
        # print str(patient) + " removed"

def load_diatance_matrix(distance_file_route):
    resi_distance_matrix_file = file(distance_file_route)
    line = resi_distance_matrix_file.readline().strip()
    resi_index = line.split("\t")[1:]
    size = len(resi_index)
    resi_distance_matrix = []
    line = resi_distance_matrix_file.readline().strip()
    while line:
        distance = []
        distances = line.split("\t")[1:]
        for dis in distances:
            distance.append(int(dis))
        resi_distance_matrix.append(distance)
        line = resi_distance_matrix_file.readline().strip()

    resi_distance_matrix_file.close()
    return resi_distance_matrix,size, resi_index

def clustering(resi_distance_matrix, size, resi_index,distance_threshold,P,R,F,cluster_truth,clusters_result_file):
    seed = random.randint(1,200)
    random.seed = seed
    cluster_set = []
    global_patient_set = set()
    order = [i for i in range(size)]
    shuffle(order)

    for i in order:
        patient = resi_index[i]
        cluster = set()
        if patient not in global_patient_set:
            cluster.add(patient)
            global_patient_set.add(patient)
            for j in range(0, size):
                distance = int(resi_distance_matrix[i][j])
                if distance < distance_threshold and resi_index[j] not in global_patient_set:
                    cluster.add(resi_index[j])
                    global_patient_set.add(resi_index[j])
            if len(cluster) > 1:
                check_convergence(resi_index, cluster, patient, global_patient_set, resi_distance_matrix)
        if len(cluster) > 0:
            cluster_set.append(cluster)

    cluster_numbers = get_cluster_count(cluster_set, cluster_truth)

    for cluster in cluster_set:
        line = ""
        for patient in cluster:
            line += str(patient) + " "
        line += "\n"
        clusters_result_file.writelines(line)
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

    TP, TN, FP_TP, total = claculate_TPTN(cluster_numbers)
    TN /= 2
    FN = combination(total) - FP_TP - TN

    # print "TP:" + str(TP)
    # print "TN:" + str(TN)
    # print "FP_TP:" + str(FP_TP)
    # print "FN:" + str(FN)
    # print total

    # RI = (TP+TN)/combination(total)
    p = TP / FP_TP
    r = TP / (TP + FN)
    f = 2 * p * r / (p + r)
    # print RI
    P.append(p)
    R.append(r)
    F.append(f)


for isresi in ["resi","unresi"]:
    resi_distance_matrix, size, resi_index = load_diatance_matrix("./output_new/distance_"+isresi+".txt")

    clusters_result_statistics_file = file("./output_new/cluster_"+isresi+"_statistcs.txt", "w")

    clusters_truth_file = file("./output_new/0.26_100_0.01_0.018_0.05_0.75_"+isresi+".txt")
    cluster_truth = get_cluster_truth(clusters_truth_file)

    clusters_result_file = file("./output_new/cluster_"+isresi+".txt", "w")

    if isresi == "resi":
        distance_threshold = [10, 30, 50, 70]
    else:
        distance_threshold = [30, 60, 90, 120]
    for threshold in distance_threshold:
        P = []
        R = []
        F = []
        for i in range(50):
            clustering(resi_distance_matrix, size, resi_index, threshold, P, R, F,cluster_truth,clusters_result_file)

        clusters_result_statistics_file.writelines("Distance_Threshlod:" + str(threshold) + "\n\n")

        line = " ".join(str(e) for e in P)
        clusters_result_statistics_file.writelines("precision: " + line + "\n")
        P_pd = pandas.Series(P)
        clusters_result_statistics_file.writelines("precision mean: " + str(P_pd.mean()) + "\n")
        clusters_result_statistics_file.writelines("precision std: " + str(P_pd.std()) + "\n\n")

        line = " ".join(str(e) for e in R)
        clusters_result_statistics_file.writelines("recall: " + line + "\n")
        R_pd = pandas.Series(R)
        clusters_result_statistics_file.writelines("recall mean: " + str(R_pd.mean()) + "\n")
        clusters_result_statistics_file.writelines("recall std: " + str(R_pd.std()) + "\n\n")

        line = " ".join(str(e) for e in F)
        clusters_result_statistics_file.writelines("F1: " + line + "\n")
        F_pd = pandas.Series(F)
        clusters_result_statistics_file.writelines("F1 mean: " + str(F_pd.mean()) + "\n")
        clusters_result_statistics_file.writelines("F1 std: " + str(F_pd.std()) + "\n\n")

    clusters_truth_file.close()
    clusters_result_file.close()
    clusters_result_statistics_file.close()
class patient:
    def __init__(self):
        self.removal = False
        self.resistant = []
        self.mutation_pos = []
        self.mutation_nuc = []
        self.resistant_mutation = set()

def load_positions(position_path):
    position_file=file(position_path)
    position_list=[]
    position=position_file.readline().strip()
    while position:
        position_list.append(position)
        position = position_file.readline().strip()
    position_file.close()
    return position_list

def load_resistant_SNPs(resi_path, SNP_positions,TB_sequence):
    resitant_SNPs_file=file(resi_path)
    line = resitant_SNPs_file.readline().strip()
    resistant_SNPs = {}
    while line:
        if line.startswith("Mycobacterium"):
            fields=line.split("\t")
            position=int(SNP_positions.index(fields[1]))
            mutation=TB_sequence[position]+" "+fields[4].upper()
            if position in resistant_SNPs.keys():
                resistant_SNPs[position] = resistant_SNPs[position]+","+mutation
            else:
                resistant_SNPs[position] = mutation
        line = resitant_SNPs_file.readline().strip()
    resitant_SNPs_file.close()
    return resistant_SNPs

def ISresistant(sequence,resistant_SNPs):
    isresi = False
    for key in resistant_SNPs.keys():
        if isresi:
            break
        Nu_seq = sequence[key]
        Mu_resis = resistant_SNPs[key].split(",")
        for Mu_resi in Mu_resis:
            Nu_resi = Mu_resi.split()[1]
            if Nu_resi == Nu_seq:
                print key
                isresi = True
    return isresi

def distance(patients,i,j):
    resistant_mutation_intesection = patients[i].resistant_mutation.intersection(patients[j].resistant_mutation)
    if len(resistant_mutation_intesection) == 0:
        return -1
    else:
        mutation_set_i = set(patients[i].mutation_pos)
        mutation_set_j = set(patients[j].mutation_pos)
        mutation_difference = mutation_set_i.symmetric_difference(mutation_set_j)
        return len(mutation_difference)

def distance_unresi(patients, i, j):
    mutation_set_i = set(patients[i].mutation_pos)
    mutation_set_j = set(patients[j].mutation_pos)
    mutation_difference = mutation_set_i.symmetric_difference(mutation_set_j)
    return len(mutation_difference)

def calculate_distance(patients,resi_index,isresi):
    number_patients = len(resi_index)
    matrix = [[0 for i in range(number_patients)] for j in range(number_patients)]
    for i in range(number_patients):
        for j in range(i+1,number_patients):
            if isresi:
                matrix[i][j] = matrix[j][i] = distance(patients,i,j)
            else:
                matrix[i][j] = matrix[j][i] = distance_unresi(patients, i, j)
    return matrix

def reconstruct_patients_list(resi_patients,unresi_patients, simulation_file):
    resi_patients_index = []
    unresi_patients_index = []
    text = simulation_file.readlines()
    i = 0
    patient_text = text[7 + 6 * i:13 + 6 * i]
    while patient_text:
        pat = patient()
        if not len(patient_text) < 6:
            mutations = patient_text[1].strip()[11:-1].split(";")
            for mutation in mutations:
                pos, nuc = mutation.split(":")
                pat.mutation_pos.append(pos)
                pat.mutation_nuc.append(nuc)
            index = int(patient_text[0].strip()[9:])
            pat.resistant = patient_text[3][14:-2]
            if len(pat.resistant) != 0:
                resi_patients.append(pat)
                resi_patients_index.append(index)
            else:
                unresi_patients.append(pat)
                unresi_patients_index.append(index)
        else:
            break

        i += 1
        pat.resistant = patient_text[3].strip()[14:-1].split()
        pat.resistant_mutation = set(patient_text[2].strip()[25:-2].split())
        patient_text = text[7 + 6 * i:13 + 6 * i]
        # pat.eventlist = eventlist_str2list(patient_text[1].strip())
        # pat.mutation = mutation_str2dic(patient_text[2].strip())
        # pat.resistant_mutation = resistant_mutation_str2set(patient_text[3].strip())
        # pat.resistant = resistant_str2bool(patient_text[4].strip())
        # pat.removal = removal_str2bool(patient_text[5].strip())
        # patients.append(pat)

    return resi_patients_index, unresi_patients_index

ancestor = file("ancestor.fasta")
TB_sequence=ancestor.readline().strip()
SNP_positions = load_positions("mutate_SNPs.txt")
resistant_SNPs = load_resistant_SNPs("resi.vcf", SNP_positions, TB_sequence)


resi_patients = []
unresi_patients = []
simulation_file = file("./output_new/simulation_unremoved0.26_100_0.01_0.0018_0.05_0.75_7.txt")
resi_index,unresi_index = reconstruct_patients_list(resi_patients,unresi_patients, simulation_file)
distance_matrix_resi = calculate_distance(resi_patients,resi_index,True)
distance_matrix_unresi = calculate_distance(unresi_patients,unresi_index,False)

resi_distance_output_file = file("./output_new/distance_resi.txt","w")
resi_index_file = file("./output_new/index_resi.txt","w")
unresi_distance_output_file = file("./output_new/distance_unresi.txt","w")
unresi_index_file = file("./output_new/index_unresi.txt","w")

index_out = "*"+"\t"+str(resi_index[0])
for i in range(1,len(resi_index)):
    index_out += "\t"+str(resi_index[i])
resi_distance_output_file.writelines(index_out+"\n")

for i in range(len(resi_index)):
    line = str(resi_index[i]) + "\t"+str(distance_matrix_resi[i][0])
    for j in range(1,len(resi_index)):
        line += "\t"+str(distance_matrix_resi[i][j])
    resi_distance_output_file.writelines(line+"\n")

for index in resi_index:
    resi_index_file.writelines(str(index)+"\n")

index_out = "*"+"\t"+str(unresi_index[0])
for i in range(1,len(unresi_index)):
    index_out += "\t"+str(unresi_index[i])
unresi_distance_output_file.writelines(index_out+"\n")

for i in range(len(unresi_index)):
    line = str(unresi_index[i]) + "\t"+str(distance_matrix_unresi[i][0])
    for j in range(1,len(unresi_index)):
        line += "\t"+str(distance_matrix_unresi[i][j])
    unresi_distance_output_file.writelines(line+"\n")

for index in unresi_index:
    unresi_index_file.writelines(str(index)+"\n")
# distance_matrix_unresi = calculate_distance(unresi_patients,unresi_index,False)
# for i in range(len(unresi_index)):
#     line = str(distance_matrix_unresi[i][0])
#     for j in range(1,len(unresi_index)):
#         line += " "+str(distance_matrix_unresi[i][j])
#     unresi_distance_output_file.writelines(line+"\n")

simulation_file.close()
resi_distance_output_file.close()
unresi_distance_output_file.close()
resi_index_file.close()



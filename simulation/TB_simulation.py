__author__ = 'Fang'
__date__= '2016.6.18'
__email__= 'fza34@sfu.ca'
__function__ = 'Tuberculosis simulation'


import random
import copy
import sys

class patient:
    def __init__(self):
        self.removal = False
        self.resistant = False
        self.eventlist = []
        self.mutation = {}
        self.resistant_mutation = set()

def encoding_nucleotide(nucleotide):
    if nucleotide == 'T':
        return 0
    elif nucleotide == 'C':
        return 1
    elif nucleotide == 'A':
        return 2
    else:
        return 3

def decoding_nucleotide(nucleotide_code):
    if nucleotide_code == 0:
        return 'T'
    elif nucleotide_code == 1:
        return 'C'
    elif nucleotide_code == 2:
        return 'A'
    else:
        return 'G'

def nucleotide_muatation(nucleotide_encoding, mutation_spectrum):
    chance=random.uniform(0,1)
    if chance <=mutation_spectrum[nucleotide_encoding][0]:
        return decoding_nucleotide(0)
    elif chance <=mutation_spectrum[nucleotide_encoding][0]+mutation_spectrum[nucleotide_encoding][1]:
        return decoding_nucleotide(1)
    elif chance <=mutation_spectrum[nucleotide_encoding][0]+mutation_spectrum[nucleotide_encoding][1]+mutation_spectrum[nucleotide_encoding][2]:
        return decoding_nucleotide(2)
    else:
        return decoding_nucleotide(3)

def load_positions(position_path):
    position_file=file(position_path)
    position_list=[]
    position=position_file.readline().strip()
    while position:
        position_list.append(position)
        position = position_file.readline().strip()
    position_file.close()
    return position_list

def get_distance_matrix(sampling_patients,TB90):
    matrix = [[0 for col in range(len(TB90))] for row in range(len(sampling_patients))]
    for i in range(0,len(sampling_patients)):
        for j in range(0,len(TB90)):
            try:
                matrix[i][j] = sum(1 for a, b in zip(sampling_patients[i], TB90[j]) if a != b)
            except:
                print str(len(sampling_patients)) + " " + str(len(TB90))
    return matrix

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

def calculate_events_number(m, n, miu, t, beta, alfa, Pt, Pr, gama):
    n_mutation = m*miu*t
    n_transmission = m*beta*t
    n_resistant = m*alfa*Pt*Pr*t
    n_removal = m*gama*t
    n_tatal = n_mutation + n_transmission + n_resistant + n_removal

    return n_tatal, n_mutation, n_transmission, n_resistant, n_removal

def if_mutate_back(mutate_position, patient, nucleotide_muatated):
    if not mutate_position in patient.mutation:
        return False
    elif not nucleotide_muatated == patient.mutation[mutate_position][0]:
        return False
    else:
        return True

def if_SNP_resistant(mutate_position, TB_sequence,  nucleotide_muatated, resistant_SNPs):
    nucleotide=TB_sequence[mutate_position]
    SNPs = resistant_SNPs[mutate_position].split(",")
    for SNP in SNPs:
        mutation_pair = SNP.split(" ")
        if nucleotide == mutation_pair[0] and nucleotide_muatated == mutation_pair[1]:
            return True
    return False

def if_mutate_resistant(mutate_position, TB_sequence,  nucleotide_muatated, resistant_SNPs):
    if not mutate_position in resistant_SNPs:
        return False
    elif not if_SNP_resistant(mutate_position, TB_sequence,  nucleotide_muatated, resistant_SNPs):
        return False
    return True

def mutationmap2string(mutation):
    mutation_str = ""
    for mutation_position in mutation.keys():
        mutation_str += str(mutation_position)+":"+mutation[mutation_position]+";"
    return mutation_str

def patients_sampling(patients, kappa, TB_sequence):
    samples = []
    unremoved = set()
    for i in range(0, len(patients)):
        if not patients[i].removal:
            unremoved.add(i)
    sample_number=int(kappa*len(unremoved))
    patients_sampling = random.sample(unremoved, sample_number)
    for i in patients_sampling:
        sequence_list = list(TB_sequence)
        for key in patients[i].mutation:
             sequence_list[key]=patients[i].mutation[key][1]
        samples.append(sequence_list)
    return samples

def mutate(patient, TB_sequence, n, mutation_spectrum,resistant_SNPs, allow_collision):
    # if not mutate_position in patient.mutation: # no mutation at this position
    #     nucleotide_origin=TB_sequence[mutate_position]
    # else: # mutation at this position
    #     nucleotide_origin=patient.mutation[mutate_position].split(" ")[1]
    if allow_collision:
        mutate_position = random.randint(0, n - 1)  # mutation
        if mutate_position in patient.mutation.keys():
            return
    else:
        mutate_positions = set(range(n)) - set(patient.mutation.keys())-set(resistant_SNPs.keys())
        if mutate_positions:
            mutate_position = random.sample(mutate_positions, 1)[0]
        else:
            return
    # print mutate_position
    nucleotide_origin = TB_sequence[mutate_position]
    nucleotide_encoding=encoding_nucleotide(nucleotide_origin)
    nucleotide_muatated=nucleotide_muatation(nucleotide_encoding, mutation_spectrum) # randomly get the mutated nucleotide according to the mutation_spectrum
    # while if_mutate_resistant(mutate_position, TB_sequence, nucleotide_muatated, resistant_SNPs): # no mutate back and no resistant mutation
        # print "TRUE"
        # print resistant_SNPs[mutate_position]
        # nucleotide_muatated = nucleotide_muatation(nucleotide_encoding, mutation_spectrum)
        # print mutate_position
    patient.mutation[mutate_position] = TB_sequence[mutate_position]+nucleotide_muatated #
    patient.eventlist.append(str(1)+" "+str(mutate_position)+" "+nucleotide_origin+" "+nucleotide_muatated) # record
    # print " "+nucleotide_origin+" "+nucleotide_muatated+"\n"

def get_resistance(patient, TB_sequence, resistant_SNPs, allow_collision):
    if allow_collision:
        mutate_position = resistant_SNPs.keys()[random.randint(0, len(resistant_SNPs) - 1)]  # mutation
        if mutate_position in patient.resistant_mutation:
            return
    else:
        mutate_positions = set(resistant_SNPs.keys()) - patient.resistant_mutation
        if mutate_positions:
            mutate_position = random.sample(mutate_positions, 1)[0]
        else:
            return
    SNPs=resistant_SNPs[mutate_position].split(",") # get the mutations at this position
    SNP_resistant=SNPs[random.randint(0, len(SNPs)-1)].split(" ")[1] # get a randomly selected mutation
    if not mutate_position in patient.mutation: # no mutation in this position
        origin_nucleotide = TB_sequence[mutate_position]
    elif not SNP_resistant == patient.mutation[mutate_position][1]: # SNP is not the mutated one
        origin_nucleotide = patient.mutation[mutate_position][1]
    # else:
    #     origin_nucleotide = "A"
    # print mutate_position
    patient.mutation[mutate_position] = TB_sequence[mutate_position]+SNP_resistant
    patient.eventlist.append(str(3)+" "+ str(mutate_position) +" "+origin_nucleotide + " " + SNP_resistant)
    patient.resistant = True
    patient.resistant_mutation.add(mutate_position)
    # print " " + origin_nucleotide + " " + SNP_resistant+"\n"

def develop_simulation(patients, TB_sequence, generation, n, mutation_spectrum, resistant_SNPs, allow_collision):
    for i in range(0, generation):
        for j in range(0, len(patients)):
            patients.append(copy.deepcopy(patients[j]))# duplication
        for j in range(0, len(patients)):
            mutate(patients[j], TB_sequence, n, mutation_spectrum, resistant_SNPs, allow_collision)


def breakout_simulation(patients, TB_sequence, n, mutation_spectrum, events_count, resistant_SNPs,allow_collision):
    P_mutation=events_count[1]/float(events_count[0])
    P_transmission=(events_count[1]+events_count[2])/float(events_count[0])
    P_resistance=(events_count[1]+events_count[2]+events_count[3])/float(events_count[0])

    for i in range(0, int(events_count[0])):
        patient=random.randint(0,len(patients)-1) # randomly select one patient
        ran = random.random() # a randon between 0-n_tatal-1
        # print patient
        if not patients[patient].removal:
            if ran <= P_mutation:
                # print "mutation at "
                mutate(patients[patient], TB_sequence, n, mutation_spectrum, resistant_SNPs, allow_collision)
            elif ran <= P_transmission:
                # print "transmit to "
                transmitted_patient=random.randint(0,len(patients)-1) # randomly select one patient to be transmitted
                while patient == transmitted_patient or patients[transmitted_patient].removal:
                    transmitted_patient=random.randint(0,len(patients)-1)
                patients[patient].eventlist.append(str(2)+" "+str(patient)+" "+str(transmitted_patient))
                patients[transmitted_patient].eventlist.append(str(2) + " " + str(patient) + " " + str(transmitted_patient)+" "+  mutationmap2string(patients[transmitted_patient].mutation))
                patients[transmitted_patient].mutation = copy.deepcopy(patients[patient].mutation)
                patients[transmitted_patient].resistant=copy.deepcopy(patients[patient].resistant)
                patients[transmitted_patient].resistant_mutation = copy.deepcopy(patients[patient].resistant_mutation)
            elif ran <= P_resistance:
                # print "resist at "
                get_resistance(patients[patient], TB_sequence, resistant_SNPs, allow_collision)
                # resistance
            else:
                patients[patient].removal = True
                patients[patient].eventlist.append(str(4))

def trace_back(patients, n, resistant_SNPs):
    patient=patients[n]
    while patient.eventlist:
        event = patient.eventlist[-1].split(" ")
        if event[0] == '4':
            patient.removal = False
            patient.eventlist.pop()
        elif event[0]=='2':
            if event[1] == str(n):
                patient.eventlist.pop()
            else:
                patient.mutation.clear()
                mutations = event[3].split(";")
                mutations = mutations[0:len(mutations)-1]
                for mutation in mutations:
                    fields=mutation.split(":")
                    patient.mutation[int(fields[0])]=fields[1]
                patient.resistant_mutation.clear()
                patient.resistant=False
                for mutate_position in  patient.mutation.keys():
                    nucleotide_muatated=list(patient.mutation[mutate_position])[1]
                    if if_mutate_resistant(mutate_position, TB_sequence, nucleotide_muatated, resistant_SNPs):
                        patient.resistant_mutation.add(mutate_position)
                        patient.resistant=True
                patient.eventlist.pop()
        elif event[0]=='1':
            SNP = list(patient.mutation[int(event[1])])
            if event[2] == SNP[0] and event[3] == SNP[1]:
                patient.mutation.pop(int(event[1]))
                patient.eventlist.pop()
            else:
                print "Bug!"
                break
        else:
            SNP = list(patient.mutation[int(event[1])])
            if event[3] == SNP[1]:
                if event[2] == SNP[0]:
                    patient.mutation.pop(int(event[1]))
                else:
                    SNP[1] = event[2]
                    patient.mutation[int(event[1])] = "".join(SNP)
                try:
                    patient.resistant_mutation.remove(int(event[1]))
                except Exception, e:
                    # print patient.eventlist[-1]
                    # print event[1]+event[2]+event[3]
                    print Exception, ":", e
                    # print if_mutate_resistant(int(event[1]), TB_sequence, event[3], resistant_SNPs)
                    # print resistant_SNPs[int(event[1])]
                if not patient.resistant_mutation:
                    patient.resistant = False
                patient.eventlist.pop()
            else:
                print "Bug!"
                break
    if patient.mutation or patient.removal or patient.resistant:
        " failed traced back!"

if __name__ == '__main__':
    TB90 = []
    TB90_file=file("./TB90_6936.txt")
    TB90_sequence=TB90_file.readline().strip()
    while TB90_sequence:
        TB90.append(list(TB90_sequence))
        TB90_sequence = TB90_file.readline().strip()
    TB90_file.close()

    generation =10  # generation in develop period
    m=pow(2,generation)# number of patients
    miu=0.26# mutation rate

    t=1 # time span 10, 15, 20, 25, 30, 35, 40, 45, 50
    beta=0 # contact/reinfection rate 0.05, 0.1, 0.15,0.2, 0.25, 0.3, . . . , 0.5
    alfa=0.3 # rate of breakdown 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5
    Pt=0.2 # probability of seeking treatment 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8
    Pr=0.3 # probability of resistant 0.01,0.05, 0.1, 0.15, 0.2, 0.3, .0.4, 0.5, 0.6, 0.7, 1
    gama=0.01 #rate of removal 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7

    kappa=0.1
    allow_collision=False # True: skip mutation event when collision False: find another position to mutate unless there are no unmutated positions


    # read ancestral sequence 6936
    # ancestor = file(sys.argv[1])
    ancestor = file("ancestor.fasta")
    TB_sequence=ancestor.readline().strip()
    ancestor.close()
    n = len(TB_sequence)  # length of TB genome sequence
    # TB_sequence="TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAGCAAAGGGCTTGGCTCAATCTCGTCCAGCCATTGACCATCGTCGAGGGGTTTGCTCTGTTATCCGTGCCGAG"
    # TB_sequence=TB_sequence[0:n]

    events_count = calculate_events_number(m, n, miu, t, beta, alfa, Pt, Pr, gama)  # calculate events numbers

    # SNP_positions = load_positions(sys.argv[2])
    SNP_positions = load_positions("mutate_SNPs.txt")
    # load drug-resistant mutations
    resistant_SNPs = load_resistant_SNPs("resi.vcf", SNP_positions, TB_sequence) # resistant SNPs
    # resistant_SNPs = {}
    # resistant_SNPs[3]="A C,A G"
    # resistant_SNPs[5]="C A"
    # for key in resistant_SNPs.keys():
    #     if key > n-1:
    #         resistant_SNPs.pop(key)

    mutation_spectrum=[[0, 0.6265, 0.1371, 0.2364], [0.5587, 0, 0.2313, 0.21], [0.1428, 0.2699, 0, 0.5873], [0.2348, 0.2272, 0.538, 0]] # mutation spectrum of nucleotides
    # mutation_spectrum = [[0, 0.33, 0.33, 0.34], [0.33, 0, 0.33, 0.34], [0.33, 0.33, 0, 0.34],[0.33, 0.33, 0.34, 0]]
    patients = [] # list of patients
    original_patient = patient()
    patients.append(original_patient)

    develop_simulation(patients, TB_sequence, generation, n, mutation_spectrum, resistant_SNPs, allow_collision)

    breakout_simulation(patients, TB_sequence, n, mutation_spectrum, events_count, resistant_SNPs, allow_collision)

    sampling_patients=patients_sampling(patients, kappa,TB_sequence)

    distance_matrix=get_distance_matrix(sampling_patients,TB90)
    # for i in range(0, m):
    #     print i
    #     trace_back(patients,i,resistant_SNPs)

    output_path = "./output/simulation_"+str(miu)+"_"+str(t)+"_"+str(beta)+"_"+str(alfa)+"_"+str(Pt)+"_"+str(Pr)+"_"+str(gama)+"_"+".txt"
    output=file(output_path, 'w')
    output.write("generation =" + str(generation))
    output.write("\tnumber of patients =" + str(m))
    output.write("\tmutation rate =" + str(miu))
    output.write("\ttime span =" + str(t))
    output.write("\tcontact/reinfection rate =" + str(beta))
    output.write("\trate of breakdown =" + str(alfa))
    output.write("\tprobability of seeking treatment =" + str(Pt))
    output.write("\tprobability of resistant =" + str(Pr))
    output.write("\trate of removal =" + str(gama))
    output.write("\trate of sampling =" + str(kappa))
    output.write("\tallow_collision =" + str(allow_collision))
    output.write("\n\nEvents number:" + str(events_count))
    output.write("\n\nPatient list: \n")
    i = 0
    for patient in patients:
        output.write("Patient: " + str(i) + "\n")
        output.write("".join(patient.eventlist)+"\n")
        output.write("mutations: "+str(patient.mutation)+"\n")
        output.write("resistant_mutations:" + str(patient.resistant_mutation)+"\n")
        output.write("isresistant: "+ str(patient.resistant)+"\n")
        output.write("isremoved: " + str(patient.removal)+"\n")
        output.write("\n")
        output.write("Distance matrix:\n")
        i += 1

    d = " \t"
    for i in range(0, len(TB90)):
        d = d + str(i) + "\t"
    d += "\n"
    output.write(d)
    for i in range(0, len(sampling_patients)):
        d = str(i) + "\t"
        for j in range(0, len(TB90)):
            d = d + str(distance_matrix[i][j]) + "\t"
        d = d +"\n"
        output.write(d)

    # print events_count[0]
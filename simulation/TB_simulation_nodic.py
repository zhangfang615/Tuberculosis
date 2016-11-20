from __future__ import division
__author__ = 'Fang'
__date__= '2016.10.21'
__email__= 'fza34@sfu.ca'
__function__ = 'Tuberculosis simulation'


import random
import copy
import numpy

class patient:
    def __init__(self):
        self.removal = False
        self.resistant = []
        self.mutation_pos = []
        self.mutation_nuc = []
        self.resistant_mutation = set()

# nucleotide to number
def encoding_nucleotide(nucleotide):
    if nucleotide == 'T':
        return 0
    elif nucleotide == 'C':
        return 1
    elif nucleotide == 'A':
        return 2
    else:
        return 3

# number to nucleotide
def decoding_nucleotide(nucleotide_code):
    if nucleotide_code == 0:
        return 'T'
    elif nucleotide_code == 1:
        return 'C'
    elif nucleotide_code == 2:
        return 'A'
    else:
        return 'G'

# based on mutation spectrum, random the mutation
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

# load mutated positions
def load_positions(position_path):
    position_file=file(position_path)
    position_list=[]
    position=position_file.readline().strip()
    while position:
        position_list.append(position)
        position = position_file.readline().strip()
    position_file.close()
    return position_list

# calculate distance matrix
def get_distance_matrix(sampling_patients,TB90):
    matrix = [[0 for col in range(len(TB90))] for row in range(len(sampling_patients))]
    for i in range(0,len(sampling_patients)):
        for j in range(0,len(TB90)):
            try:
                matrix[i][j] = sum(1 for a, b in zip(sampling_patients[i], TB90[j]) if a != b)
            except:
                print str(len(sampling_patients)) + " " + str(len(TB90))
    return matrix

# load resistant SNPs
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

# calculate event numbers
def calculate_events_number(m, n, miu, t, beta, P_resi, gama):
    n_mutation = m*miu*t
    n_transmission = m*beta*t
    n_resistant = m*P_resi*t
    n_removal = m*gama*t
    n_tatal = n_mutation + n_transmission + n_resistant + n_removal
    return n_tatal, n_mutation, n_transmission, n_resistant, n_removal

# judge whether a SNP is resistant SNP
def if_SNP_resistant(mutate_position, TB_sequence,  nucleotide_muatated, resistant_SNPs):
    nucleotide=TB_sequence[mutate_position]
    SNPs = resistant_SNPs[mutate_position].split(",")
    for SNP in SNPs:
        mutation_pair = SNP.split(" ")
        if nucleotide == mutation_pair[0] and nucleotide_muatated == mutation_pair[1]:
            return True
    return False

# judge whether a mutation is resistant
def if_mutate_resistant(mutate_position, TB_sequence,  nucleotide_muatated, resistant_SNPs):
    if not mutate_position in resistant_SNPs:
        return False
    elif not if_SNP_resistant(mutate_position, TB_sequence,  nucleotide_muatated, resistant_SNPs):
        return False
    return True


def mutationmap2string(mutation_pos,mutation_nuc):
    mutation_str = ""
    for i in range(0,len(mutation_pos)):
        mutation_str += str(mutation_pos[i]) + ":"+mutation_nuc[i]+";"
    return mutation_str

# sampling un removed patients
def patients_sampling(patients, kappa, TB_sequence):
    samples = []
    unremoved = set()
    for i in range(0, len(patients)):
        if not patients[i].removal:
            unremoved.add(i)
    # print str(len(unremoved))
    sample_number=int(kappa*len(unremoved))
    patients_sampling = random.sample(unremoved, sample_number)
    return patients_sampling


def patients_sample(unremoved, kappa):
    sample_number=int(kappa*len(unremoved))
    patients_sampling = random.sample(unremoved, sample_number)
    return patients_sampling

# mutation
def mutate(patient, TB_sequence, n, mutation_spectrum,resistant_SNPs, allow_collision):
    if allow_collision:
        mutate_position = random.randint(0, n - 1)  # mutation
        if mutate_position in patient.mutation_pos:
            return
    else:
        mutate_positions = set(range(n)) - set(patient.mutation_pos)-set(resistant_SNPs.keys())
        if mutate_positions:
            mutate_position = random.sample(mutate_positions, 1)[0]
        else:
            return
    nucleotide_origin = TB_sequence[mutate_position]
    nucleotide_encoding=encoding_nucleotide(nucleotide_origin)
    nucleotide_muatated=nucleotide_muatation(nucleotide_encoding, mutation_spectrum) # randomly get the mutated nucleotide according to the mutation_spectrum
    patient.mutation_pos.append(mutate_position)
    patient.mutation_nuc.append(TB_sequence[mutate_position]+nucleotide_muatated) #
    # patient.eventlist.append(str(1)+" "+str(mutate_position)+" "+nucleotide_origin+" "+nucleotide_muatated) # record

# resistance aquisition
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
    if not mutate_position in patient.mutation_pos: # no mutation in this position
        origin_nucleotide = TB_sequence[mutate_position]
    elif not SNP_resistant == patient.mutation_nuc[patient.mutation_pos.index(mutate_position)][1]: # SNP is not the mutated one
        origin_nucleotide = patient.mutation[mutate_position][1]

    patient.mutation_pos.append(mutate_position)
    patient.mutation_nuc.append(TB_sequence[mutate_position]+SNP_resistant)
    # patient.eventlist.append(str(3)+" "+ str(mutate_position) +" "+origin_nucleotide + " " + SNP_resistant)
    patient.resistant.append(3)
    patient.resistant_mutation.add(mutate_position)

# develop stage
def develop_simulation(patients, TB_sequence, generation, n, mutation_spectrum, resistant_SNPs, allow_collision):
    for i in range(0, generation):
        print i
        for j in range(0, len(patients)):
            patients.append(copy.deepcopy(patients[j]))# duplication
        for j in range(0, len(patients)):
            span = random.random()
            if span <=0.7:
                a = random.randrange(5,16)
            else:
                a = random.randrange(16,35)
            for count in range(0,a):
                mutate(patients[j], TB_sequence, n, mutation_spectrum, resistant_SNPs, allow_collision)

# break out stage
def breakout_simulation_sample(patients, TB_sequence, n, mutation_spectrum, events_count, resistant_SNPs,allow_collision,substitution,kappa,globalevent_list):
    P_mutation=events_count[1]/float(events_count[0])
    P_transmission=(events_count[1]+events_count[2])/float(events_count[0])
    P_resistance=(events_count[1]+events_count[2]+events_count[3])/float(events_count[0])
    unremoved = set(range(0, len(patients)))
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
                if len(patients[transmitted_patient].resistant) == 0 or len(patients[patient].resistant) > 0:
                    if len(patients[patient].resistant) == 0:
                        globalevent_list.append(str(1) + " " + str(patient) + " " + str(transmitted_patient))
                    else:
                        globalevent_list.append(str(2) + " " + str(patient) + " " + str(transmitted_patient))
                    # patients[patient].eventlist.append(str(2) + " " + str(patient) + " " + str(transmitted_patient))
                    # patients[transmitted_patient].eventlist.append(str(2) + " " + str(patient) + " " + str(transmitted_patient) + " " + mutationmap2string(patients[transmitted_patient].mutation))
                    patients[transmitted_patient].mutation_pos = copy.deepcopy(patients[patient].mutation_pos)
                    patients[transmitted_patient].mutation_nuc = copy.deepcopy(patients[patient].mutation_nuc)
                    if len(patients[patient].resistant) > 0:
                        patients[transmitted_patient].resistant.append(2)
                    patients[transmitted_patient].resistant_mutation = copy.deepcopy(patients[patient].resistant_mutation)
            elif ran <= P_resistance:
                # print "resist at "
                get_resistance(patients[patient], TB_sequence, resistant_SNPs, allow_collision)
                globalevent_list.append(str(4) + " "  + str(patient))
                # resistance
            else:
                sub = random.random()
                if sub <= substitution:
                    transmit_patient = random.randint(0, len(patients) - 1)
                    while patient == transmit_patient or patients[transmit_patient].removal or len(patients[transmit_patient].resistant)> 0:
                        transmit_patient = random.randint(0, len(patients) - 1)
                    globalevent_list.append(str(3) + " " + str(transmit_patient) + " " + str(patient))
                    # patients[transmit_patient].eventlist.append(str(2) + " " + str(transmit_patient) + " " + str(patient))
                    # patients[patient].eventlist.append(str(2) + " " + str(transmit_patient) + " " + str(patient) + " " + mutationmap2string(patients[patient].mutation))
                    patients[patient].mutation_pos = copy.deepcopy(patients[transmit_patient].mutation_pos)
                    patients[patient].mutation_nuc = copy.deepcopy(patients[transmit_patient].mutation_nuc)
                    patients[patient].resistant = []
                    patients[patient].resistant_mutation = copy.deepcopy(patients[transmit_patient].resistant_mutation)
                    patients[patient].removal = False
                else:
                    globalevent_list.append(str(5) + " " + str(patient))   
                    patients[patient].removal = True
                    # patients[patient].eventlist.append(str(4))

def TB_output(patients,m, miu, t, beta, P_resi,gama,output_path):
    output = file(output_path, 'w')
    output.write("generation =" + str(generation))
    output.write("\tnumber of patients =" + str(m))
    output.write("\tmutation rate =" + str(miu))
    output.write("\ttime span =" + str(t))
    output.write("\tcontact/reinfection rate =" + str(beta))
    output.write("\tprobability of resistant =" + str(P_resi))
    output.write("\trate of removal =" + str(gama))
    output.write("\trate of sampling =" + str(kappa))
    output.write("\tallow_collision =" + str(allow_collision))
    output.write("\n\nEvents number:" + str(events_count))
    output.write("\n\nRandom seed:" + str(seed))
    output.write("\n\nPatient list: \n")
    count_resistant = 0
    count_unremoved = 0
    count_unremoved_resistant = 0
    i = 0
    for pat in patients:
        output.write("Patient: " + str(i) + "\n")
        # output.write(eventlist2string(pat.eventlist).strip() + "\n")
        output.write("mutations: " + mutationmap2string(pat.mutation_pos,pat.mutation_nuc)+ "\n")
        output.write("resistant_mutations:" + str(pat.resistant_mutation) + "\n")
        output.write("isresistant: " + str(pat.resistant) + "\n")
        output.write("isremoved: " + str(pat.removal) + "\n")
        if len(pat.resistant) > 0:
            count_resistant += 1
        if not pat.removal:
            count_unremoved += 1
        if len(pat.resistant) > 0 and not pat.removal:
            count_unremoved_resistant += 1
        output.write("\n")
        i += 1
    output.write("\nresistant patients number" + str(count_resistant))
    output.write("\nunremoved patients number" + str(count_unremoved))
    output.write("\nunremoved resistant patients number" + str(count_unremoved_resistant))
    output.close()

def final_sequences_output(patients, unremoved, output_sequences_path,TB_sequence):
    output_sequences = file(output_sequences_path,"w")
    for patient in unremoved:
        sequence = list(TB_sequence)
        for i in range(0,len(patients[patient].mutation_pos)):
            mutated_nuc = patients[patient].mutation_nuc[i][-1]
            sequence[patients[patient].mutation_pos[i]] = mutated_nuc
        sequence_output = "".join(a for a in sequence)+"\n"
        output_sequences.writelines(sequence_output)
    output_sequences.close()

def TB_unremoved_output(patients,m, miu, t, beta, P_resi,gama,output_path, unremoved_path):
    output = file(output_path, 'w')
    output.write("generation =" + str(generation))
    output.write("\tnumber of patients =" + str(m))
    output.write("\tmutation rate =" + str(miu))
    output.write("\ttime span =" + str(t))
    output.write("\tcontact/reinfection rate =" + str(beta))
    output.write("\tprobability of resistant =" + str(P_resi))
    output.write("\trate of removal =" + str(gama))
    output.write("\trate of sampling =" + str(kappa))
    output.write("\tallow_collision =" + str(allow_collision))
    output.write("\n\nEvents number:" + str(events_count))
    output.write("\n\nRandom seed:" + str(seed))
    output.write("\n\nPatient list: \n")

    unremoved_output = file(unremoved_path, 'w')
    unremoved_output.write("generation =" + str(generation))
    unremoved_output.write("\tnumber of patients =" + str(m))
    unremoved_output.write("\tmutation rate =" + str(miu))
    unremoved_output.write("\ttime span =" + str(t))
    unremoved_output.write("\tcontact/reinfection rate =" + str(beta))
    unremoved_output.write("\tprobability of resistant =" + str(P_resi))
    unremoved_output.write("\trate of removal =" + str(gama))
    unremoved_output.write("\trate of sampling =" + str(kappa))
    unremoved_output.write("\tallow_collision =" + str(allow_collision))
    unremoved_output.write("\n\nEvents number:" + str(events_count))
    unremoved_output.write("\n\nRandom seed:" + str(seed))
    unremoved_output.write("\n\nPatient list: \n")

    unremoved = set()

    count_resistant = 0
    count_unremoved = 0
    count_unremoved_resistant = 0
    i = 0
    for pat in patients:
        output.write("Patient: " + str(i) + "\n")
        # output.write(eventlist2string(pat.eventlist).strip() + "\n")
        output.write("mutations: " + mutationmap2string(pat.mutation_pos,pat.mutation_nuc)+ "\n")
        output.write("resistant_mutations:" + str(pat.resistant_mutation) + "\n")
        output.write("isresistant: " + str(pat.resistant) + "\n")
        output.write("isremoved: " + str(pat.removal) + "\n")
        if len(pat.resistant) > 0:
            count_resistant += 1
        if not pat.removal:
            unremoved_output.write("Patient: " + str(i) + "\n")
            unremoved_output.write("mutations: " + mutationmap2string(pat.mutation_pos, pat.mutation_nuc) + "\n")
            unremoved_output.write("resistant_mutations:" + str(pat.resistant_mutation) + "\n")
            unremoved_output.write("isresistant: " + str(pat.resistant) + "\n")
            unremoved_output.write("isremoved: " + str(pat.removal) + "\n")
            unremoved_output.write("\n")
            unremoved.add(i)
            count_unremoved += 1
        if len(pat.resistant) > 0 and not pat.removal:
            count_unremoved_resistant += 1
        output.write("\n")
        i += 1
    output.write("\nresistant patients number" + str(count_resistant))
    output.write("\nunremoved patients number" + str(count_unremoved))
    output.write("\nunremoved resistant patients number" + str(count_unremoved_resistant))
    unremoved_output.write("\nunremoved patients number" + str(count_unremoved))
    unremoved_output.write("\nunremoved resistant patients number" + str(count_unremoved_resistant))
    output.close()
    unremoved_output.close()
    return unremoved

def TB_statistics(patients,unremoved,kappa,statistics_path):
    statistics_output = file(statistics_path,'w')
    for Ka in kappa:
        count_sampling = []
        count_resistant = []
        count_unresistant = []
        count_resistant_trans = []
        count_resistant_acq = []
        count_resistant_mevent = []
        for i in range(0,5):
            seed = random.randint(0, 100000)
            random.Random(seed)
            sampling_patients = patients_sample(unremoved, Ka)
            count_resistant.append(0)
            count_unresistant.append(0)
            count_resistant_trans.append(0)
            count_resistant_acq.append(0)
            count_resistant_mevent.append(0)

            count_sampling.append(len(sampling_patients))
            sampling_patients.sort()
            resistant_sample = set()

            for sample in sampling_patients:  # count resistant samples and separate resistant and unresistant samples
                if patients[sample].resistant:
                    count_resistant[i] += 1
                    resistant_sample.add(sample)
                else:
                    count_unresistant[i] += 1

            for sample in resistant_sample:
                if patients[sample].resistant[-1] == 2:
                    count_resistant_trans[i] += 1
                else:
                    count_resistant_acq[i] += 1

                if len(patients[sample].resistant) > 1:
                    count_resistant_mevent[i] += 1
        count_sampling_mean = numpy.mean(count_sampling)

        count_resistant_mean = numpy.mean(count_resistant)
        count_resistant_std = numpy.std(count_resistant)
        resistant_ratio = count_resistant_mean / count_sampling_mean * 100

        count_resistant_trans_mean = numpy.mean(count_resistant_trans)
        count_resistant_trans_std = numpy.std(count_resistant_trans)
        resistant_trans_ratio = count_resistant_trans_mean / count_sampling_mean * 100

        count_resistant_acq_mean = numpy.mean(count_resistant_acq)
        count_resistant_acq_std = numpy.std(count_resistant_acq)
        resistant_acq_ratio = count_resistant_acq_mean / count_sampling_mean * 100

        count_resistant_mevent_mean = numpy.mean(count_resistant_mevent)
        count_resistant_mevent_std = numpy.std(count_resistant_mevent)
        if not count_resistant_mean == 0:
            resistant_mevent_ratio = count_resistant_mevent_mean / count_resistant_mean * 100
        else:
            resistant_mevent_ratio = 0

        output_line =     "Sampling rate: "+str(Ka) \
                          + "\n" + str(count_sampling_mean) \
                          + "\t" + str(count_resistant_mean) + "\t" + str(count_resistant_std) + "\t" + str(resistant_ratio) \
                          + "\t" + str(count_resistant_trans_mean) + "\t" + str(count_resistant_trans_std) + "\t" + str(resistant_trans_ratio) \
                          + "\t" + str(count_resistant_acq_mean) + "\t" + str(count_resistant_acq_std) + "\t" + str(resistant_acq_ratio) \
                          + "\t" + str(count_resistant_mevent_mean) + "\t" + str(count_resistant_mevent_std) + "\t" + str(resistant_mevent_ratio) \
                          + "\n"
        statistics_output.writelines(output_line)
    statistics_output.close()

if __name__ == '__main__':
    TB90 = []
    TB90_file=file("./TB90_6936.txt")
    TB90_sequence=TB90_file.readline().strip()
    while TB90_sequence:
        TB90.append(list(TB90_sequence))
        TB90_sequence = TB90_file.readline().strip()
    TB90_file.close()

    generation =20 # generation in develop period
    m=pow(2,generation) # number of patients
    miu=0.26 # mutation rate

    t_list=[100] # time span 10, 15, 20, 25, 30, 35, 40, 45, 50
    beta_list = [0.01] # transmission rate
    P_resi_list = [0.004] # acquisition rate
    gama_list = [0.05] # removal rate
    substitution_list = [0.75] # substitution rate for removed patients
    kappa=[0.001,0.005,0.01] # sampling rate
    allow_collision=False # True: skip mutation event when collision False: find another position to mutate unless there are no unmutated positions

    ancestor = file("ancestor.fasta")
    TB_sequence=ancestor.readline().strip()
    ancestor.close()
    n = len(TB_sequence)  # length of TB genome sequence
    # TB_sequence="TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAGCAAAGGGCTTGGCTCAATCTCGTCCAGCCATTGACCATCGTCGAGGGGTTTGCTCTGTTATCCGTGCCGAG"
    # TB_sequence=TB_sequence[0:n]

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

    # mutation_spectrum=[[0, 0.6265, 0.1371, 0.2364], [0.5587, 0, 0.2313, 0.21], [0.1428, 0.2699, 0, 0.5873], [0.2348, 0.2272, 0.538, 0]] # mutation spectrum of nucleotides
    mutation_spectrum = [[0, 0.33, 0.33, 0.34], [0.33, 0, 0.33, 0.34], [0.33, 0.33, 0, 0.34],[0.33, 0.33, 0.34, 0]]
    globalevent_list = []
    for t in t_list:
        for beta in beta_list:
            for P_resi in P_resi_list:
                for gama in gama_list:
                    for substitution in substitution_list:
                        output_path = "./output_new/simulation_" + str(miu) + "_" + str(t) + "_" + str(beta) + "_" + str(P_resi) + "_" + str(gama) + "_" + str(substitution)
                        output_unremoved_path = "./output_new/simulation_unremoved" + str(miu) + "_" + str(t) + "_" + str(beta) + "_" + str(P_resi) + "_" + str(gama) + "_" + str(substitution)
                        output_eventlsit_path = "./output_new/eventlist" + str(miu) + "_" + str(t) + "_" + str(beta) + "_" + str(P_resi) + "_" + str(gama) + "_" + str(substitution)+".txt"
                        output_sequences_path = "./output_new/final_sequences.txt"
                        if True:
                            seed = random.randint(0, 100000)
                            random.Random(seed)
                            patients = []  # list of patients
                            original_patient = patient()
                            patients.append(original_patient)
                            develop_simulation(patients, TB_sequence, generation, n, mutation_spectrum, resistant_SNPs,allow_collision) # Develop stage

                            for sample_time in range(1,3):
                                t, P_resi, gama= 20,0,0.001
                                events_count = calculate_events_number(m, n, miu, t, beta, P_resi,gama)  # calculate events numbers
                                # first 40 years in break stage with no resistance aquisitions (trans rate = 0.01, P_resi = 0, removal rate gama = 0.001)
                                breakout_simulation_sample(patients, TB_sequence, n, mutation_spectrum,events_count, resistant_SNPs, allow_collision,substitution, kappa,globalevent_list)
                                O_path = output_path + "_"+str(sample_time)+".txt"
                                # TB_output(patients,m,  miu, t, beta, P_resi,gama,O_path)

                            # breakout_simulation(patients, TB_sequence, n, mutation_spectrum, events_count,resistant_SNPs, allow_collision, substitution)
                            for sample_time in range(3,8):
                                t, P_resi,gama = 20, 0.004,0.05
                                events_count = calculate_events_number(m, n, miu, t, beta, P_resi, gama)
                                # 100 years oin break stage (trans rate = 0.01, P_resi = 0.004, removal rate gama = 0.05)
                                breakout_simulation_sample(patients, TB_sequence, n, mutation_spectrum, events_count,resistant_SNPs, allow_collision, substitution, kappa,globalevent_list)
                                O_path = output_path + "_" + str(sample_time) + ".txt"
                                unremoved_path = output_unremoved_path + "_" + str(sample_time) + ".txt"
                                statistics_path = output_unremoved_path + "_" + str(sample_time) + "_statistics.txt"
                                # unremoved = TB_unremoved_output(patients, m,  miu, t, beta, P_resi, gama, O_path,unremoved_path)
                                # TB_statistics(patients,unremoved,kappa,statistics_path)

                            unremoved = TB_unremoved_output(patients, m, miu, t, beta, P_resi, gama, O_path,unremoved_path) # outpout
                            TB_statistics(patients, unremoved, kappa, statistics_path) # statistics

                            # output evnetlist
                            eventlist_file = file(output_eventlsit_path,"w")
                            for event in globalevent_list:
                                eventlist_file.writelines(event+"\n")
                            eventlist_file.close()

                            # outout final sequences
                            final_sequences_output(patients, unremoved, output_sequences_path,TB_sequence)
                            del patients
                            del original_patient

__author__ = 'Fang'
__date__= '2016.6.18'
__email__= 'fza34@sfu.ca'
__function__ = 'Tuberculosis simulation'


import random
import copy

class patient:
    def __init__(self):
        self.removal = False
        self.resistant = False
        self.eventlist = []
        self.mutation = {}
        self.resistant_mutation = set()

def encoding_nucleotide(nucleotide):
    if nucleotide == 'A':
        return 0
    elif nucleotide == 'T':
        return 1
    elif nucleotide == 'C':
        return 2
    else:
        return 3

def decoding_nucleotide(nucleotide_code):
    if nucleotide_code == 0:
        return 'A'
    elif nucleotide_code == 1:
        return 'T'
    elif nucleotide_code == 2:
        return 'C'
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
        print "".join(sequence_list)

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
        mutate_positions = set(range(n)) - set(patient.mutation.keys())
        if mutate_positions:
            mutate_position = random.sample(mutate_positions, 1)[0]
        else:
            return
    nucleotide_origin = TB_sequence[mutate_position]
    nucleotide_encoding=encoding_nucleotide(nucleotide_origin)
    nucleotide_muatated=nucleotide_muatation(nucleotide_encoding, mutation_spectrum) # randomly get the mutated nucleotide according to the mutation_spectrum
    while if_mutate_resistant(mutate_position, TB_sequence, nucleotide_muatated, resistant_SNPs): # no mutate back and no resistant mutation
        nucleotide_muatated = nucleotide_muatation(nucleotide_encoding, mutation_spectrum)
    patient.mutation[mutate_position] = TB_sequence[mutate_position]+nucleotide_muatated #
    patient.eventlist.append(str(1)+" "+str(mutate_position)+" "+nucleotide_origin+" "+nucleotide_muatated) # record

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
    patient.mutation[mutate_position] = TB_sequence[mutate_position]+SNP_resistant
    patient.eventlist.append(str(3)+" "+ str(mutate_position) +" "+origin_nucleotide + " " + SNP_resistant)
    patient.resistant = True
    patient.resistant_mutation.add(mutate_position)

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
        if not patients[patient].removal:
            if ran <= P_mutation:
                mutate(patients[patient], TB_sequence, n, mutation_spectrum, resistant_SNPs, allow_collision)
            elif ran <= P_transmission:
                transmitted_patient=random.randint(0,len(patients)-1) # randomly select one patient to be transmitted
                while patient == transmitted_patient or patients[transmitted_patient].removal:
                    transmitted_patient=random.randint(0,len(patients)-1)
                patients[patient].eventlist.append(str(2)+" "+str(patient)+" "+str(transmitted_patient))
                patients[transmitted_patient].eventlist.append(str(2) + " " + str(patient) + " " + str(transmitted_patient)+" "+  mutationmap2string(patients[transmitted_patient].mutation))
                patients[transmitted_patient].mutation = copy.deepcopy(patients[patient].mutation)
                patients[transmitted_patient].resistant=copy.deepcopy(patients[patient].resistant)
                patients[transmitted_patient].resistant_mutation = copy.deepcopy(patients[patient].resistant_mutation)
            elif ran <= P_resistance:
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
                    print Exception, ":", e
                if not patient.resistant_mutation:
                    patient.resistant = False
                patient.eventlist.pop()
            else:
                print "Bug!"
                break
    if patient.mutation or patient.removal or patient.resistant:
        " failed traced back!"

if __name__ == '__main__':
    generation = 10  # generation in develop period
    m=pow(2,generation)# number of patients
    n=200# length of TB genome sequence
    miu=0.26 # mutation rate
    t=50 # time span
    beta=0.05 # contact/reinfection rate
    alfa=0.3 # rate of breakdown
    Pt=0.5 # probability of seeking treatment
    Pr=0.3 # probability of resistant
    gama=0.01 #rate of removal
    kappa=0.1
    allow_collision=False # True: skip mutation event when collision False: find another position to mutate unless there are no unmutated positions

    events_count = calculate_events_number(m, n, miu, t, beta, alfa, Pt, Pr, gama) # calculate events numbers

    TB_sequence="TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAGCAAAGGGCTTGGCTCAATCTCGTCCAGCCATTGACCATCGTCGAGGGGTTTGCTCTGTTATCCGTGCCGAG"
    TB_sequence=TB_sequence[0:n]
    mutation_spectrum=[[0, 0.33, 0.33, 0.34], [0.33, 0, 0.33, 0.34], [0.33, 0.33, 0, 0.34], [0.33, 0.33, 0.34, 0]] # mutation spectrum of nucleotides

    resistant_SNPs = {} # resistant SNPs
    resistant_SNPs[3]="A C,A G"
    resistant_SNPs[5]="C A"
    for key in resistant_SNPs.keys():
        if key > n-1:
            resistant_SNPs.pop(key)

    patients = [] # list of patients
    original_patient = patient()
    patients.append(original_patient)

    develop_simulation(patients, TB_sequence, generation, n, mutation_spectrum, resistant_SNPs, allow_collision)

    breakout_simulation(patients, TB_sequence, n, mutation_spectrum, events_count, resistant_SNPs, allow_collision)

    patients_sampling(patients, kappa,TB_sequence)

    for i in range(0, m):
        trace_back(patients,i,resistant_SNPs)
    print events_count[0]
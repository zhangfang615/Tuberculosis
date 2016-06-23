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
    n_mutation = m*n*miu*t
    n_transmission = m*beta*t
    n_resistant = m*alfa*Pt*Pr*t
    n_removal = m*gama*t
    n_tatal = n_mutation + n_transmission + n_resistant + n_removal

    return n_tatal, n_mutation, n_transmission, n_resistant, n_removal

def if_mutate_back(mutate_position, patient, nucleotide_muatated):
    if not mutate_position in patient.mutation:
        return False
    elif not nucleotide_muatated == patient.mutation[mutate_position].split(" ")[0]:
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

def mutate(patient, TB_sequence, mutation_spectrum, mutate_position, resistant_SNPs):
    if not mutate_position in patient.mutation: # no mutation at this position
        nucleotide_origin=TB_sequence[mutate_position]
    else: # mutation at this position
        nucleotide_origin=patient.mutation[mutate_position].split(" ")[1]
    nucleotide_encoding=encoding_nucleotide(nucleotide_origin)
    nucleotide_muatated=nucleotide_muatation(nucleotide_encoding, mutation_spectrum) # randomly get the mutated nucleotide according to the mutation_spectrum
    if not if_mutate_back(mutate_position, patient, nucleotide_muatated) and not if_mutate_resistant(mutate_position, TB_sequence, nucleotide_muatated, resistant_SNPs): # no mutate back and no resistant mutation
        patient.mutation[mutate_position] = TB_sequence[mutate_position]+" "+nucleotide_muatated #
        patient.eventlist.append(str(1)+" "+str(mutate_position)+" "+nucleotide_origin+" "+nucleotide_muatated) # record

def get_resistance(patient, TB_sequence, resistant_SNPs):
    mutate_position = resistant_SNPs.keys()[random.randint(0, len(resistant_SNPs)-1)] # get a randomly selected mutate position
    SNPs=resistant_SNPs[mutate_position].split(",") # get the mutations at this position
    SNP_resistant=SNPs[random.randint(0, len(SNPs)-1)].split(" ")[1] # get a randomly selected mutation
    if not mutate_position in patient.mutation: # no mutation in this position
        origin_nucleotide = TB_sequence[mutate_position]
        patient.mutation[mutate_position] = origin_nucleotide+" "+SNP_resistant
        patient.eventlist.append(str(3)+" "+ str(mutate_position) +" "+origin_nucleotide + " " + SNP_resistant)
        patient.resistant = True
    elif not SNP_resistant == patient.mutation[mutate_position].split(" ")[1]: # SNP is not the mutated one
        origin_nucleotide = patient.mutation[mutate_position].split(" ")[0]
        patient.mutation[mutate_position] =origin_nucleotide +" "+SNP_resistant
        patient.eventlist.append(str(3)+" "+ str(mutate_position) +" "+origin_nucleotide + " " + SNP_resistant)
        patient.resistant = True

def develop_simulation(patients, TB_sequence, generation, n, mutation_spectrum, resistant_SNPs):
    for i in range(0, generation):
        for j in range(0, len(patients)):
            patients.append(copy.deepcopy(patients[j]))# duplication
        for j in range(0, len(patients)):
            mutate_position=random.randint(0,n-1) # mutation
            mutate(patients[j], TB_sequence, mutation_spectrum, mutate_position, resistant_SNPs)


def breakout_simulation(patients, TB_sequence, n, mutation_spectrum, events_count, resistant_SNPs):
    P_mutation=events_count[1]/float(events_count[0])
    P_transmission=(events_count[1]+events_count[2])/float(events_count[0])
    P_resistance=(events_count[1]+events_count[2]+events_count[3])/float(events_count[0])

    for i in range(0, int(events_count[0])):
        patient=random.randint(0,len(patients)-1) # randomly select one patient
        ran = random.random() # a randon between 0-n_tatal-1
        if not patients[patient].removal:
            if ran <= P_mutation:
                mutate_position=random.randint(0,n-1) # mutation
                mutate(patients[patient], TB_sequence, mutation_spectrum, mutate_position, resistant_SNPs)
            elif ran <= P_transmission:
                transmitted_patient=random.randint(0,len(patients)-1) # randomly select one patient to be transmitted
                while patient == transmitted_patient:
                    transmitted_patient=random.randint(0,len(patients)-1)
                if not patients[transmitted_patient].removal:
                    patients[patient].eventlist.append(str(2)+" "+str(patient)+" "+str(transmitted_patient))
                    patients[transmitted_patient].mutation = patients[patient].mutation
                    patients[transmitted_patient].eventlist.append(str(2)+" "+str(patient)+" "+str(transmitted_patient))
            elif ran <= P_resistance:
                get_resistance(patients[patient], TB_sequence, resistant_SNPs)
                # resistance
            else:
                patients[patient].removal = True
                patients[patient].eventlist.append(str(4))

if __name__ == '__main__':
    m=32# number of patients
    n=200# length of TB genome sequence
    miu=0.0005 # mutation rate
    t=50 # time span
    beta=0.2 # contact/reinfection rate
    alfa=0.3 # rate of breakdown
    Pt=0.5 # probability of seeking treatment
    Pr=0.3 # probability of resistant
    gama=0.01 #rate of removal

    generation=5 # generation in develop period
    events_count = calculate_events_number(m, n, miu, t, beta, alfa, Pt, Pr, gama) # calculate events numbers

    TB_seqeunce="TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTTAACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAGCAAAGGGCTTGGCTCAATCTCGTCCAGCCATTGACCATCGTCGAGGGGTTTGCTCTGTTATCCGTGCCGAG"

    mutation_spectrum=[[0, 0.33, 0.33, 0.34], [0.33, 0, 0.33, 0.34], [0.33, 0.33, 0, 0.34], [0.33, 0.33, 0.34, 0]] # mutation spectrum of nucleotides

    resistant_SNPs = {} # resistant SNPs
    resistant_SNPs[3]="A C,A G"
    resistant_SNPs[5]="C A"

    patients = [] # list of patients
    original_patient = patient()
    patients.append(original_patient)

    develop_simulation(patients, TB_seqeunce, generation, n, mutation_spectrum, resistant_SNPs)

    breakout_simulation(patients, TB_seqeunce, n, mutation_spectrum, events_count, resistant_SNPs)

    print events_count[0]

from __future__ import division
__author__ = 'Fang'
__date__= '2016.9.16'
__email__= 'fza34@sfu.ca'
__function__ = 'Tuberculosis simulation statistics'

import os
import random
import ast

class patient:
    def _init_(self):
        self.removal = False
        self.resistant = False
        self.eventlist = []
        self.mutation = {}
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

def eventlist_str2list(eventlist_string):
    eventlist = []
    events = eventlist_string.split("\t")
    for event in events:
        eventlist.append(event)
    return eventlist

def mutation_str2dic(mutation_string):
    mutation_string = mutation_string[11:]
    return ast.literal_eval(mutation_string)

def  resistant_mutation_str2set(resistant_mutation_string):
    resistant_mutation_string = resistant_mutation_string[20:]
    # resistant_mutation_string= "set([2352, 1094, 1066, 3659, 3632, 3708])"
    return eval(resistant_mutation_string)

def resistant_str2bool(resistant_string):
    resistant_string = resistant_string.split(' ')[1]
    if resistant_string == "True":
        return True
    else:
        return False

def removal_str2bool(removal_string):
    removal_string = removal_string.split(' ')[1]
    if removal_string == "True":
        return True
    else:
        return False

def reconstruct_patients_list(patients, simulation_file):
    text = simulation_file.readlines()
    for i in range(1024):
        pat = patient()
        patient_text=text[7+7*i:14+7*i]
        pat.eventlist = eventlist_str2list(patient_text[1].strip())
        pat.mutation = mutation_str2dic(patient_text[2].strip())
        pat.resistant_mutation = resistant_mutation_str2set(patient_text[3].strip())
        pat.resistant = resistant_str2bool(patient_text[4].strip())
        pat.removal = removal_str2bool(patient_text[5].strip())
        patients.append(pat)
    return patients

def patients_sampling(unremoved, kappa):
    sample_number=int(kappa*len(unremoved))
    patients_sampling = random.sample(unremoved, sample_number)
    return patients_sampling

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

def get_resistant_eventlist(patients,n, resistant_SNPs):
    resistant_eventlist = []
    patient = patients[n]
    while patient.eventlist:
        event = patient.eventlist[-1].split(" ")
        if event[0] == '4':
            patient.removal = False
            patient.eventlist.pop()
        elif event[0] == '2':
            if event[1] == str(n):
                patient.eventlist.pop()
            else:
                patient.mutation.clear()
                mutations = event[3].split(";")
                mutations = mutations[0:len(mutations) - 1]
                for mutation in mutations:
                    fields = mutation.split(":")
                    patient.mutation[int(fields[0])] = fields[1]
                patient.resistant_mutation.clear()
                res = patient.resistant # record resistant
                patient.resistant = False
                for mutate_position in patient.mutation.keys():
                    nucleotide_muatated = list(patient.mutation[mutate_position])[1]
                    if if_mutate_resistant(mutate_position, TB_sequence, nucleotide_muatated, resistant_SNPs):
                        patient.resistant_mutation.add(mutate_position)
                        patient.resistant = True
                if res == True and patient.resistant ==False:
                    resistant_eventlist.append('2')
                patient.eventlist.pop()
        elif event[0] == '1':
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
                resistant_eventlist.append('3')
                if not patient.resistant_mutation:
                    patient.resistant = False
                patient.eventlist.pop()
            else:
                print "Bug!"
                break
    if patient.mutation or patient.removal or patient.resistant:
        " failed traced back!"
    return resistant_eventlist

if __name__ == '__main__':
    for filename in os.listdir('E:/PYTHON_PROJECTS/TB_simulation/test'):
        ancestor = file("E:/PYTHON_PROJECTS/TB_simulation/ancestor.fasta")
        TB_sequence = ancestor.readline().strip()
        ancestor.close()
        SNP_positions = load_positions("E:/PYTHON_PROJECTS/TB_simulation/mutate_SNPs.txt")
        resistant_SNPs = load_resistant_SNPs("E:/PYTHON_PROJECTS/TB_simulation/resi.vcf", SNP_positions, TB_sequence)

        kappa = 0.1
        patients = []
        parameters = filename[0:len(filename)-4].split('_')[2:6]
        simulation_file = file("E:/PYTHON_PROJECTS/TB_simulation/output/"+filename)
        reconstruct_patients_list(patients, simulation_file)

        unremoved = set()
        for i in range(0, len(patients)):
            if not patients[i].removal:
                unremoved.add(i)

        count_resistant = []
        count_resistant_trans = []
        count_resistant_acq = []
        count_resistant_mevent = []
        count_unresistant_oncere = []

        for i in range (0,10):
            seed = random.randint(0, 100000)
            random.Random(seed)
            sampling_patients = patients_sampling(unremoved, kappa)

            sampling_patients.sort()
            resistant_sample=set()
            unresistant_sample=set()

            for sample in sampling_patients: # count resistant samples and separate resistant and unresistant samples
                if patients[sample].resistant:
                    count_resistant[i] += 1
                    resistant_sample.add(sample)
                else:
                    unresistant_sample.add(sample)

            for sample in resistant_sample:
                resistant_eventlist = get_resistant_eventlist(patients, sample, resistant_SNPs)
                if resistant_eventlist[0] == '2':
                    count_resistant_trans[i] += 1
                else:
                    count_resistant_acq[i] += 1

                if len(resistant_eventlist) > 1:
                    count_resistant_mevent[i] += 1

            for sample in unresistant_sample:
                resistant_eventlist = get_resistant_eventlist(patients, sample, resistant_SNPs)
                if len(resistant_eventlist) >0:
                    count_unresistant_oncere[i] += 1
        print 0
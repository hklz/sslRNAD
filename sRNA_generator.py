from nupack import *
import random
import re
import sys
import subprocess
import argparse
import primer3


class webInput():
    def __init__(self, Number, scaffold, score, NED, DTD, FsRNA, Structure, MFE):
        self.Number = Number
        self.Core_scaffold = scaffold
        self.Score = score
        self.Normalized_Ensemble_Defect = NED
        self.Detail_Thermodynamics = DTD
        self.Final_sRNA_Sequence = FsRNA
        self.Final_sRNA_Structure = Structure
        self.MFE = MFE

def random_nucleotide_generator():
    bases = ['A','U','G','C']
    return bases[random.randint(0,3)]

def reverse_nucleotide(nt):
    if nt == 'T':
        nt = 'U'
    bases = 'AUGC'
    reverse_bases = 'UACG'
    return reverse_bases[bases.find(nt)]
#print(reverse_nucleotide('A'))

def seq_rev(seq):
    temp = []
    rev_seq = ''
    for i in seq:
        temp.append(i)
    for j in range(0,len(temp)):
        rev_seq = rev_seq + reverse_nucleotide(temp[len(temp)-1-j])
    return rev_seq

def Seqrev3_5(seq):
    Seq = seq_rev(seq)
    seq_list = list(Seq)
    temp = ''
    for i in seq_list:
        temp = temp + i
    return temp

def random_sRNA_seq_generator():
    sRNA_seq =''
    #Randomly generate Forward seq
    for i in range(0,11):
        sRNA_seq = sRNA_seq + random_nucleotide_generator()
    stem = sRNA_seq[0:7]
    sRNA_seq = sRNA_seq+seq_rev(stem)
    #print(sRNA_seq)
    return sRNA_seq

def sRNA_design_ensemble_defect_evl(sRNAseq, target_structure = ): #degenerate nucleotide codes N: A,C,G,U
    sRNA_1 = Domain(sRNAseq, name='sRNA_1')
    sRNA = TargetStrand([sRNA_1], name='sRNA')
    Complex_7stem_4loop = TargetComplex([sRNA], target_structure, name='')
    t1 = TargetTube(on_targets={Complex_7stem_4loop: 1e-8}, name = 't1')
    my_model = Model()
    my_tubes = [t1]
    my_design = tube_design(tubes=my_tubes,     hard_constraints=[], soft_constraints=[],  defect_weights=None, options=None, model=my_model)
    my_results = my_design.run(trials=1)[0].defects.ensemble_defect
    return my_results

def RNA_eval_analysis(Seq, target_structure = ''):
    file = open('RNA_eval_input.txt', 'w')
    file.write(Seq + '\n' + target_structure)
    file.close()
    ret = subprocess.Popen(["RNAeval", "-v", 'RNA_eval_input.txt'], stdout=subprocess.PIPE)
    stdout_var = ret.communicate()
    data_content = stdout_var[0].decode('utf-8')
    print(data_content)
    pattern = re.compile(r, re.S)
    data_loop_thermo = pattern.findall(data_content)
    MFE = float((data_content[len(data_content)-8:len(data_content)-2]).replace(" ", ""))
    print(MFE)
    #print(data_loop_thermo)acronyms
    detail_thermo = []
    for i in range(1, len(data_loop_thermo)):
        #print(int(data_loop_thermo[i]))
        detail_thermo.append(int(data_loop_thermo[i]))
    return detail_thermo

def RNA_fold_final_structure_analysis(Seq = ''):
    file = open('RNA_fold_input.fa', 'w')
    file.write(Seq + '\n')
    file.close()
    ret = subprocess.Popen(["RNAfold", "-i", 'RNA_fold_input.fa'], stdout=subprocess.PIPE)
    stdout_var = ret.communicate()
    data_content = stdout_var[0].decode('utf-8')
    #print(data_content[len(Seq) + 1:len(Seq) + 1 + len(Seq)])
    #return (data_content[len(Seq)+1:len(Seq)+1+len(Seq)])
    return (data_content[len(Seq)+1:])


def generater_desire_seq(Target_sequence = ,Repression_level = 'S', trial = 5, defect = 0.07): #Target_sequence = 'ATGCAGTCATCGTAGCAGTCAGTC'
    counter = 0
    Target_sequence = Target_sequence.upper()
    #Text_result = ""
    webInput_result = []
    while (counter<trial):
        random_test_sRNA = random_sRNA_seq_generator()
        #print(random_test_sRNA)
        Seq_defect = sRNA_design_ensemble_defect_evl(random_test_sRNA)
        if Seq_defect <= defect:
            final_structrue = RNA_fold_final_structure_analysis(seq_rev(Target_sequence)+random_test_sRNA+'UUUUUUUUU')
            #print(final_structrue)
            #print(final_structrue[len(Target_sequence):len(Target_sequence)+27])
            if final_structrue[len(Target_sequence):len(Target_sequence)+27] == '':
                detail_thermo = RNA_eval_analysis(random_test_sRNA)
                predition_result = predition_func(detail_thermo)
                if Repression_level == 'S' and round(predition_result)/() >=60:
                    counter = counter + 1
                    # Target_sequence + random_test_sRNA+'UUUUUUUUU'
                    Score = str(round((predition_result -)
                    DTD = '['
                    for i in detail_thermo[0:7]:
                        DTD = DTD + str(i) + ', '
                    DTD = DTD[:len(DTD) - 2]
                    DTD = DTD + ']'
                    MFE = (((final_structrue[len(Target_sequence)+28:]).replace('\n',' ')).replace('(','')).replace(')','')+'kcal/mol'
                    result_input = webInput(counter, random_test_sRNA, Score, Seq_defect, DTD, seq_rev(Target_sequence)+random_test_sRNA+'UUUUUUUUU', final_structrue[:len(Target_sequence)+27], MFE)
                    DTD = '['
                    webInput_result.append(result_input.__dict__)
                    webInput_json_dict = {"webInputdata": webInput_result}

                elif Repression_level == 'M' and round(predition_resuland round(predition_result:
                    counter = counter + 1
                    Score = str(round((predition_result - )
                    DTD = '['
                    for i in detail_thermo[0:7]:
                        DTD = DTD + str(i) + ', '
                    DTD = DTD[:len(DTD) - 2]
                    DTD = DTD + ']'
                    MFE = (((final_structrue[len(Target_sequence)+28:]).replace('\n',' ')).replace('(','')).replace(')','')+'kcal/mol'
                    result_input = webInput(counter, random_test_sRNA, Score, Seq_defect, DTD, seq_rev(Target_sequence)+random_test_sRNA+'UUUUUUUUU', final_structrue[:len(Target_sequence)+27], MFE)
                    DTD = '['
                    webInput_result.append(result_input.__dict__)
                    webInput_json_dict = {"webInputdata": webInput_result}

                elif Repression_level == 'W' and round(predition_result0) <20:
                    counter = counter + 1
                    Score = str(round((predition_result - ))
                    DTD = '['
                    for i in detail_thermo[0:7]:
                        DTD = DTD + str(i) + ', '
                    DTD = DTD[:len(DTD) - 2]
                    DTD = DTD + ']'
                    MFE = (((final_structrue[len(Target_sequence)+28:]).replace('\n',' ')).replace('(','')).replace(')','')+'kcal/mol'
                    result_input = webInput(counter, random_test_sRNA, Score, Seq_defect, DTD, seq_rev(Target_sequence)+random_test_sRNA+'UUUUUUUUU', final_structrue[:len(Target_sequence)+27], MFE)
                    DTD = '['
                    webInput_result.append(result_input.__dict__)
                    webInput_json_dict = {"webInputdata": webInput_result}
    return webInput_json_dict



print(generater_desire_seq())

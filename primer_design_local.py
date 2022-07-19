from Bio import SeqIO
import json
import methods.sRNA_generator as TsR
import methods.primer_design_test as PrD

class fileResult():
    def __init__(self, name, Final_sRNA_sequence, PrimerF_1, PrimerF_2, PrimerR_1, PrimerR_2):
        self.name = "Anti-%s sRNA" % name
        self.Final_sRNA_sequence = Final_sRNA_sequence
        self.PrimerF_1 = PrimerF_1
        self.PrimerF_2 = PrimerF_2
        self.PrimerR_1 = PrimerR_1
        self.PrimerR_2 = PrimerR_2

def targets_read(path = "test.fasta"):
    target_dict = {}
    for seq_record in SeqIO.parse(path, "fasta"):
        #print(seq_record.id)
        #print(type(seq_record.seq))
        target_dict[seq_record.id] = str(seq_record.seq)
    return target_dict

def inputFromfile(path, left_arm, right_arm): #promoter_default as J23119
    targetsFromfile = targets_read(path)
    valueTotarget = {v:k for k,v in targetsFromfile.items()}
    #read file
    """
    result = ""
    for i in targetsFromfile.values():
        result = result + valueTotarget[i] + '\n' + TsR.generater_desire_seq_file_input(Target_sequence = i) + "\n"
    print(result)
    """
    result = []
    for i in targetsFromfile.values():
        primer_design =  PrD.primer_design(left_arm, i, right_arm)
        #string output
        #result = result + ("Anti-%s sRNA:\nFinal sRNA sequence: %s\nAuto_design Primers:\n%s") % (valueTotarget[i], sRNA_seq, primer_design) + "\n"
        #object array output
        single_output = fileResult(valueTotarget[i], i, primer_design[0],primer_design[1],primer_design[2], primer_design[3])
        result.append(single_output.__dict__)
    fileInput_json_dict = {"FileInputdata":result}
    return fileInput_json_dict

print(inputFromfile("erg_sRNA.fasta", "CGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATGCTGGATCC", "AAGCTTAGATCTATTACCCTGTTATCCCTACTCGAGTTCATGTGCAGCTCCATCAGCA"))
import primer3
from Bio import SeqIO
class fileResult():
    def __init__(self, name, Final_sRNA_sequence, PrimerF_1, PrimerF_2, PrimerR_1, PrimerR_2):
        self.name = "Anti-%s sRNA" % name
        self.Final_sRNA_sequence = Final_sRNA_sequence
        self.PrimerF_1 = PrimerF_1
        self.PrimerF_2 = PrimerF_2
        self.PrimerR_1 = PrimerR_1
        self.PrimerR_2 = PrimerR_2

def seq_rev(targetSequence):
    rev_seq =""
    Seq = targetSequence[::-1].upper()
    for i in Seq:
        if i == 'A':
            rev_seq += 'T'
        elif i == 'T':
            rev_seq += 'A'
        elif i == 'C':
            rev_seq += 'G'
        elif i == 'G':
            rev_seq += 'C'
    return rev_seq

def Tmlength(seq):
    length = 0
    # revese_seq = seq_rev(seq)
    for i in range(len(seq)):
        if primer3.calcTm(seq[:i]) >= 53 and i >=24 and i<35:
            length = i
            break
    return seq[:length]

def primer_design(left_arm, sRNA_sequence, right_arm):
    fragment = (left_arm+sRNA_sequence+right_arm).upper().replace("U" ,"T")
    fragment_rev = seq_rev(fragment)

    #sRNA_insert = sRNA_sequence.replace("U", "T")
    primer_R_1_End = Tmlength(seq_rev(left_arm).replace("U", "T"))
    primer_R_1 = seq_rev(sRNA_sequence[:45-len(primer_R_1_End)]) + primer_R_1_End
    primer_R_1 = primer_R_1.upper().replace("U", "T")

    primer_F_1_End = Tmlength(right_arm)
    primer_F_1 = sRNA_sequence[len(sRNA_sequence)+len(primer_F_1_End)-45:len(sRNA_sequence)] + primer_F_1_End
    primer_F_1 = primer_F_1.upper().replace("U", "T")


    #insert fragment 5-3
    primer_R_2_End = Tmlength(primer_R_1)
    primer_F_2_End = Tmlength(primer_F_1)

    R_position = len(fragment) - fragment_rev.find(primer_R_2_End)+1
    F_position = fragment.find(primer_F_2_End)

    mid = int((R_position + F_position)/2)

    primer_F_2 = fragment[mid-9:F_position] + primer_F_2_End
    primer_R_2 = seq_rev(fragment[R_position-1:mid+9]) + primer_R_2_End

    primer_output = [primer_F_1, primer_F_2, primer_R_1, primer_R_2]
    return primer_output

def targets_read(path = "test.fasta"):
    target_dict = {}
    for seq_record in SeqIO.parse(path, "fasta"):
        #print(seq_record.id)
        #print(type(seq_record.seq))
        target_dict[seq_record.id] = str(seq_record.seq)
    return target_dict

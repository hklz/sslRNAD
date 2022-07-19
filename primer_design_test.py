import primer3

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
    primer_R_1 = seq_rev(sRNA_sequence[:60-len(primer_R_1_End)]) + primer_R_1_End
    primer_R_1 = primer_R_1.upper().replace("U", "T")

    primer_F_1_End = Tmlength(right_arm)
    primer_F_1 = sRNA_sequence[len(sRNA_sequence)+len(primer_F_1_End)-60:len(sRNA_sequence)] + primer_F_1_End
    primer_F_1 = primer_F_1.upper().replace("U", "T")


    #insert fragment 5-3
    primer_R_2_End = Tmlength(primer_R_1)
    primer_F_2_End = Tmlength(primer_F_1)

    R_position = len(fragment) - fragment_rev.find(primer_R_2_End)+1
    F_position = fragment.find(primer_F_2_End)

    mid = int((R_position + F_position)/2)

    primer_F_2 = fragment[mid-9:F_position] + primer_F_2_End
    primer_R_2 = seq_rev(fragment[R_position-1:mid+9]) + primer_R_2_End

    #Final output
    #string output
    #primer_output = "Primer_F_1: %s\nPrimer_F_2: %s\nPrimer_R_1: %s\nPrimer_R_2: %s\n" % (primer_F_1, primer_F_2, primer_R_1, primer_R_2)
    #arryoutput
    primer_output = [primer_F_1, primer_F_2, primer_R_1, primer_R_2]
    return primer_output

    #return [primer_F_1, primer_F_2, primer_R_1, primer_R_2]
"""
left_arm = "CTTACGCATCTGTGCGGTATTTCACACCGCATATGCTGGATCC"
sRNA = "ttgacagctagctcagtcctaggtataatgctagcGGTAAGTTCCATTGGTTCAAACATGGGCGCGTCAGCGCGCCCTTTTTTTTT"
right_arm = "AAGCTTAGATCTATTACCCTGTTATCCCTACTCGAGTTCA"

print(primer_design(left_arm,sRNA,right_arm))
"""
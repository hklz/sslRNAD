from Bio import SeqIO
import json
from openpyxl import Workbook
import methods.sRNA_generator as TsR
import methods.primer_design_test as PrD
import sys, getopt

class fileResult():
    def __init__(self, name, Final_sRNA_sequence, PrimerF_1, PrimerF_2, PrimerR_1, PrimerR_2):
        self.name = "Anti-%s sRNA" % name
        self.Final_sRNA_sequence = Final_sRNA_sequence
        self.PrimerF_1 = PrimerF_1
        self.PrimerF_2 = PrimerF_2
        self.PrimerR_1 = PrimerR_1
        self.PrimerR_2 = PrimerR_2

    def __display__(self):
        print(self.name)
        print(self.Final_sRNA_sequence)
        print(self.name+"-PrimerF-1", end=":")
        print(self.PrimerF_1)
        print(self.name+"-PrimerF-2", end=":")
        print(self.PrimerF_2)
        print(self.name+"-PrimerR-1", end=":")
        print(self.PrimerR_1)
        print(self.name+"-PrimerR-2", end=":")
        print(self.PrimerR_2)


def targets_read(path = "test.fasta"):
    target_dict = {}
    for seq_record in SeqIO.parse(path, "fasta"):
        #print(seq_record.id)
        #print(type(seq_record.seq))
        target_dict[seq_record.id] = str(seq_record.seq)
    return target_dict

def flank_sequence_read(path):
    flanks = []
    for seq_record in SeqIO.parse(path, "fasta"):
        flanks.append(str(seq_record.seq))
    return flanks

def inputFromfile(path, left_arm, right_arm, promoter = "ttgacagctagctcagtcctaggtataatgctagc"): #promoter_default as J23119
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
        sRNA_seq = TsR.generater_desire_seq_file_input(Target_sequence = i[:24])
        Insert_sequence = promoter + sRNA_seq
        primer_design =  PrD.primer_design(left_arm, Insert_sequence, right_arm)
        single_output = fileResult(valueTotarget[i], sRNA_seq, primer_design[0],primer_design[1],primer_design[2], primer_design[3])
        #single_output.__display__()
        result.append(single_output.__dict__)
        print("Anti-%s sRNA is designed. Related primers are generated." % valueTotarget[i])
    fileInput_json_dict = {"FileInputdata":result}
    return fileInput_json_dict


input_flank_sequence_path = ''
input_target_gene_path = ''
output_file_name = 'sRNA_output'
Useage = "Usage: Demo_batch_design.py -s <flanks_sedquens_file> -i <target_genes_file> -o <output_file_name>"
try:
    opts, args = getopt.getopt(sys.argv[1:], "hs:i:o:", ["help", "flanks_seq_file=", "target_genes_file=", "output_file_name"])
except getopt.GetoptError:
    print(Useage)
    sys.exit(5)
for o, a in opts:
    if o == "-h":
        print(Useage)
        sys.exit(2)
    elif o == '-s':
        input_flank_sequence_path = a
    elif o == '-i':
        input_target_gene_path = a
    elif o == '-o':
        output_file_name = a


flanks = flank_sequence_read(input_flank_sequence_path)

result_dict_list = inputFromfile(input_target_gene_path, flanks[0], flanks[1], flanks[2])

rows = []
for i in result_dict_list["FileInputdata"]:
    everyRow1 = [i['name'], i['Final_sRNA_sequence'], i['name'] + "-PrimerF_1", i['PrimerF_1']] #[name, sRNA_sequence, PrimerName, Sequence]
    everyRow2 = [i['name'], i['Final_sRNA_sequence'], i['name'] + "-PrimerF_2", i['PrimerF_2']]
    everyRow3 = [i['name'], i['Final_sRNA_sequence'], i['name'] + "-PrimerR_1", i['PrimerR_1']]
    everyRow4 = [i['name'], i['Final_sRNA_sequence'], i['name'] + "-PrimerR_2", i['PrimerR_2']]
    rows.append(everyRow1)
    rows.append(everyRow2)
    rows.append(everyRow3)
    rows.append(everyRow4)

output_book = Workbook()
output_sheet = output_book.active

for i in rows:
    output_sheet.append(i)
output_book.save(output_file_name+".xlsx")
print("sRNA design is finished!")
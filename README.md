# Tierney-Python2-FinalProject
My final project portfolio for the second part of the Python course at Louisiana Tech University 

**Sequence Objects Pt 1 **

from Bio.Seq import Seq
from Bio.Seq import Seq
y_seq = Seq ("GATCG")
my_seq = Seq ("GATCG")
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
0 G
1 A
2 T
3 C
4 G
# we can also print the length of each sequence 
print(len(my_seq))
5
print(my_seq[0])
G
print(my_seq[4])
G
print(my_seq[2])
T
Seq("AAAA").count("AA")
2
my_seq = Seq ("GATCGATGGGCCTATATAGGATCGAAAATCGC")
len(my_seq)
32
my_seq.count("G")
9
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
46.875
from Bio.SeqUtils import gc_fraction
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
gc_fraction(my_seq)
0.46875
my_seq[4:12]
Seq('GATGGGCC')
my_seq[0::3]
Seq('GCTGTAGTAAG')
my_seq[1::3]
Seq('AGGCATGCATC')
my_seq[2:3]
Seq('T')
my_seq[::-1]
Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCTAG')
str(my_seq)
'GATCGATGGGCCTATATAGGATCGAAAATCGC'
fasta_format_string = ">Name\n%s\n" % my_seq
print(fasta_format_string)
>Name
GATCGATGGGCCTATATAGGATCGAAAATCGC

seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
seq1 + seq2 
Seq('ACGTAACCGG')
seq2 + seq1
Seq('AACCGGACGT')
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTTGCA")]
spacer = Seq("N" *10)
spacer.join(contigs)
Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTTGCA')
dna_seq = Seq("acgtACGT")
​
dna_seq
Seq('acgtACGT')
dna_seq.upper()
Seq('ACGTACGT')
dna_seq.lower()
Seq('acgtacgt')
dna_seq.upper()
Seq('ACGTACGT')
"gtac" in dna_seq
False
"GTAC" in dna_seq
False
dna_seq
Seq('acgtACGT')
dna_seq = dna_seq.upper()
"GTAC" in dna_seq
True
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
my_seq.complement()
Seq('CTAGCTACCCGGATATATCCTAGCTTTTAGCG')
my_seq.reverse_complement()
Seq('GCGATTTTCGATCCTATATAGGCCCATCGATC')
protein_seq = Seq("EVRNAK")
protein_seq.complement()
Seq('EBYNTM')
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
​
coding_dna
Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
template_dna = coding_dna.reverse_complement()
template_dna 
Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')
messenger_rna = coding_dna.transcribe()
messenger_rna
Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')
template_dna.reverse_complement().transcribe()
Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')
messenger_rna.back_transcribe()
Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
messenger_rna.translate()
Seq('MAIVMGR*KGAR*')
coding_dna.translate(table="Vertebrate Mitochondrial")
Seq('MAIVMGRWKGAR*')
coding_dna.translate(table=2)
Seq('MAIVMGRWKGAR*')
coding_dna.translate(to_stop = True)
Seq('MAIVMGR')
coding_dna.translate(table =2, to_stop =True)
Seq('MAIVMGRWKGAR')
coding_dna.translate(table = 2, stop_symbol = "!")
Seq('MAIVMGRWKGAR!')
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
gene.translate(table = "Bacterial")
Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HR*')
gene.translate(table = "Bacterial", to_stop = True)
Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')
gene.translate(table = "Bacterial", cds = True)
Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')
from Bio.Data import CodonTable
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
print(standard_table)
Table 1 Standard, SGC0

  |  T      |  C      |  A      |  G      |
--+---------+---------+---------+---------+--
T | TTT F   | TCT S   | TAT Y   | TGT C   | T
T | TTC F   | TCC S   | TAC Y   | TGC C   | C
T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
--+---------+---------+---------+---------+--
C | CTT L   | CCT P   | CAT H   | CGT R   | T
C | CTC L   | CCC P   | CAC H   | CGC R   | C
C | CTA L   | CCA P   | CAA Q   | CGA R   | A
C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
--+---------+---------+---------+---------+--
A | ATT I   | ACT T   | AAT N   | AGT S   | T
A | ATC I   | ACC T   | AAC N   | AGC S   | C
A | ATA I   | ACA T   | AAA K   | AGA R   | A
A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
--+---------+---------+---------+---------+--
G | GTT V   | GCT A   | GAT D   | GGT G   | T
G | GTC V   | GCC A   | GAC D   | GGC G   | C
G | GTA V   | GCA A   | GAA E   | GGA G   | A
G | GTG V   | GCG A   | GAG E   | GGG G   | G
--+---------+---------+---------+---------+--
print(mito_table)
Table 2 Vertebrate Mitochondrial, SGC1

  |  T      |  C      |  A      |  G      |
--+---------+---------+---------+---------+--
T | TTT F   | TCT S   | TAT Y   | TGT C   | T
T | TTC F   | TCC S   | TAC Y   | TGC C   | C
T | TTA L   | TCA S   | TAA Stop| TGA W   | A
T | TTG L   | TCG S   | TAG Stop| TGG W   | G
--+---------+---------+---------+---------+--
C | CTT L   | CCT P   | CAT H   | CGT R   | T
C | CTC L   | CCC P   | CAC H   | CGC R   | C
C | CTA L   | CCA P   | CAA Q   | CGA R   | A
C | CTG L   | CCG P   | CAG Q   | CGG R   | G
--+---------+---------+---------+---------+--
A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
--+---------+---------+---------+---------+--
G | GTT V   | GCT A   | GAT D   | GGT G   | T
G | GTC V   | GCC A   | GAC D   | GGC G   | C
G | GTA V   | GCA A   | GAA E   | GGA G   | A
G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
--+---------+---------+---------+---------+--
mito_table.stop_codons
['TAA', 'TAG', 'AGA', 'AGG']
mito_table.start_codons
['ATT', 'ATC', 'ATA', 'ATG', 'GTG']
seq = Seq("ACGT")
"ACGT" == seq1
True
seq1 == "ACGT"
True
unknown_seq = Seq(None, 10)
unknown_seq
Seq(None, length=10)
len(unknown_seq)
10
seq= Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
seq[1000:1020]
Seq(None, length=20)
seq[117512690:117512700]
Seq('CCTGAATGTG')
seq[117512670:]
Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)
seq=Seq ("ACGT")
undefined_seq = Seq(None, length = 10)
seq + undefined_seq + seq
Seq({0: 'ACGT', 14: 'ACGT'}, length=18)
my_seq = Seq ("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
from Bio.Seq import MutableSeq
mutable_seq = MutableSeq(my_seq)
mutable_seq
MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')
mutable_seq[5]="C"
 mutable_seq
MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')
mutable_seq.remove("T")
mutable_seq
MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')
mutable_seq.reverse()
mutable_seq
MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')
new_seq = Seq(mutable_seq)
new_seq
Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate 
my_string = "GCTGAAGCTGAATCAG"
reverse_complement(my_string)
'CTGATTCAGCTTCAGC'
back_transcribe(my_string)
'GCTGAAGCTGAATCAG'
translate(my_string)
'AEAES'


****In the Seq Objects section we learned how to use DNA, RNA, or protein sequences and how we can code to get complements, the translation, how to edit them, and reverse these sequences. 

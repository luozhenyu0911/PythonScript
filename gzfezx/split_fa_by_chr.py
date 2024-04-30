from Bio import SeqIO
import sys

if len(sys.argv)!= 2 :
    print("Usage: python split_fa_by_chr.py ref_genome.fasta")
    sys.exit(1)
    
for record in SeqIO.parse(sys.argv[1], 'fasta') :
    with open( record.id+".fasta", "a") as output_handle :
        SeqIO.write(record, output_handle, 'fasta')
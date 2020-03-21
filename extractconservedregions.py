from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

def get_residue_number_from_column_index(aln, id, col):
    rec = next((r for r in aln if r.id == id), None)
    j = 0
    rn=0
    for i, res in enumerate(rec.seq):
        if j==col:
            print("%d %s" %(j, res))
            if res=='-':
                return rn
            else:
                return rn+1
        else:
            j+=1
            if res!='-':
                rn+=1

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

#feature type could be gene, CDS, protein, UTR and so on - please refer to GenBank nomenclature
'''
def get_location_start_end_as_list(loc):
    return [int(s) for s in loc.split() if s.isdigit()]
'''
def main():
    if len(sys.argv) < 3:
        print("Usage: python3 ExtractFeatureSequenceFromGenbankGenome.py <input.gb> <featuretype> <outdir>")
        sys.exit(0)
    else:
        outdir=sys.argv[3]
        if not outdir.endswith("/"):
            outdir=outdir+"/"
        for srec in SeqIO.parse(sys.argv[1],'gb'):
            print(srec.id)
            for f in srec.features:
                sid=""
                gname=""
                proteinid=""
                if f.type == sys.argv[2]:
                    for q1 in f.qualifiers:
                        if q1 == 'gene':
                            sid=srec.id+"__"+f.type+"__"+str(f.location.start)+"__"+str(f.location.end)+"__"+f.qualifiers[q1][0]
                            gname=f.qualifiers[q1][0]
                            featrecord=SeqRecord(f.extract(srec.seq), id=sid, name="", description="")
                            fo=outdir+sid+".fa"
                            SeqIO.write(featrecord, fo, "fasta")
                        if q1 == 'protein_id':
                            proteinid=f.qualifiers[q1][0]
                        if f.type == "CDS" and q1 == 'translation':
                            proid=srec.id+"__Protein__"+proteinid+"__"+gname
                            sqq=Seq(f.qualifiers[q1][0], generic_protein,)
                            featrecord=SeqRecord(sqq, id=proid, name="", description="")
                            fo=outdir+proid+".fa"
                            SeqIO.write(featrecord, fo, "fasta")
                #print(f.extract(srec.seq))
            #print(srec.features[0])
#python3 ExtractGeneNucleotideFromGenbankGenome.py <genbankgenomesfile> gene <ouputdir>
if __name__ == "__main__":
    main()

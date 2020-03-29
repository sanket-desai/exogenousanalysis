import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Seq
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
                if f.type == sys.argv[2]:
                    print(f.location.end)
                    for q1 in f.qualifiers:
                        if q1 == 'gene':
                            sid=srec.id+"__"+f.type+"__"+str(f.location.start)+"__"+str(f.location.end)+"__"+f.qualifiers[q1][0]
                featrecord=SeqRecord(f.extract(srec.seq), id=sid, name="", description="")
                fo=outdir+sid+".fa"
                SeqIO.write(featrecord, fo, "fasta")
                #print(f.extract(srec.seq))
            #print(srec.features[0])

if __name__ == "__main__":
    main()

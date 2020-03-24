from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import sys

class ConservedCluster(object):
    def __init__(self, name, start, end, sequence):
        self.name_=name
        self.msastartindex_=start
        self.msaendindex_=end
        self.seq_=sequence
    def __str__(self):
        return "%s\t%d\t%d\t%d\t%s" %(self.name_, self.msastartindex_,self.msaendindex_, self.msaendindex_-self.msastartindex_+1 ,  self.seq_)
    def __len__(self):
        return self.msaendindex_-self.msastartindex_+1
class DNAMSAFeatures(object):
    def __init__(self, fname, format="fasta"):
        self.alignment_= AlignIO.read(open(fname),format)
        self.purines_="AG"
        self.pyrimidines_="CT"
        self.conserved_residue_list_=self.get_conserved_nucleotides_as_list()
    def get_column_as_list(self,index):
        return list(self.alignment_[:,index])
    def get_unique_elements_in_column_as_list(self,index):
        return list(set(list(self.alignment_[:,index])))
    def is_column_conserved(self,index):
        li=self.get_unique_elements_in_column_as_list(index)
        if li == None:
            print("Bad alignment file!!")
            sys.exit(0)
        if len(li)==1:
            return True
        elif len(li)==2:
            return (li[0] in self.purines_ and li[1] in self.purines_) or (li[0] in self.pyrimidines_ and li[1] in self.pyrimidines_)
        else:
            return False
    def get_column_element(self,index):
        li=self.get_unique_elements_in_column_as_list(index)
        if len(li)==1:
            return li[0]
        elif len(li)==2:
            if li[0].upper() in self.purines_ and li[1].upper() in self.purines_:
                return "R"
            elif li[0].upper() in self.pyrimidines_ and li[1].upper() in self.pyrimidines_:
                return "Y"
            else:
                return "-"
        else:
            return "-"
    def get_conserved_nucleotides_as_list(self):
        conres=[]
        cres=""
        for i in range(0, self.alignment_.get_alignment_length()):
            conres.append(self.get_column_element(i))
        return conres
    def get_sequence_number_from_column_index(self, id, col):
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
    def get_conserved_cluster_list(self):
        cluststart=-1
        clustend=-1
        cindex=0
        consclustlist=[]
        #print(self.get_conserved_nucleotides_as_list())
        for i in self.get_conserved_nucleotides_as_list():
            if i=="-" and cluststart==-1:
                pass
            elif i!="-" and cluststart==-1:
                cluststart=cindex
            elif i=="-" and cluststart>-1 and clustend==-1:
                clustend=cindex-1
                if clustend-cluststart+1 >= 1:
                    cc=ConservedCluster("clust_"+str(cluststart)+"_"+str(clustend), cluststart, clustend, "".join(self.conserved_residue_list_[cluststart:clustend+1]) )
                    print(str(cc))
                    consclustlist.append(cc)
                    cluststart=-1
                    clustend=-1
                else:
                    cluststart=-1
                    clustend=-1
            cindex+=1
        return consclustlist
    def write_conserved_clusters_to_file(self, fname):
        fo=open(fname , 'w')
        for c in self.get_conserved_cluster_list():
            fo.write(str(c))
            fo.write("\n")
        fo.close()

def main():
    if len(sys.argv)!=3:
        print("usage: DNAMSAFeatures.py msafile outputfile")
        sys.exit(0)
    else:
        d=DNAMSAFeatures(sys.argv[1])
        d.write_conserved_clusters_to_file(sys.argv[2])
        print("Compute success!!")
if __name__=="__main__":
    main()

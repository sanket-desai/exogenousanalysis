from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from difflib import SequenceMatcher
from Bio import SeqIO
import sys

class LongestCommonFragment(object):
    def __init__(self, ifa,format='fasta'):
        self.farecords_=[]
        for r in SeqIO.parse(ifa,format):
            self.farecords_.append(r)
    def lcs(self, S,T):
        m = len(S)
        n = len(T)
        counter = [[0]*(n+1) for x in range(m+1)]
        longest = 0
        lcs_set = set()
        for i in range(m):
            for j in range(n):
                if S[i] == T[j]:
                    c = counter[i][j] + 1
                    counter[i+1][j+1] = c
                    if c > longest:
                        lcs_set = set()
                        longest = c
                        lcs_set.add(S[i-c+1:i+1])
                    elif c == longest:
                        lcs_set.add(S[i-c+1:i+1])
        return lcs_set
    def lcs_seqrecords(sr1, sr2):
        naame=sr1.name+"__"+sr2.name
        if len(name) > 20:
            naame="__"+sr1[:7]
        lset=self.lcs( str(sr1.seq), str(sr2.seq))
        lsrecs=[]
        for ls in lset:
            lsrecs.append( SeqRecord( Seq(ls), id=naame, name="", description="" ) )
        return lsrecs
    #x is the min length
    def topcs(self, S,T,x):
        m = len(S)
        n = len(T)
        counter = [[0]*(n+1) for x in range(m+1)]
        longest = 0
        lcs_set = set()
        for i in range(m):
            for j in range(n):
                if S[i] == T[j]:
                    c = counter[i][j] + 1
                    counter[i+1][j+1] = c
                    if c >=x:
                        #lcs_set = set()
                        longest = c
                        lcs_set.add(S[i-c+1:i+1])
                    #elif c == longest:
                    #lcs_set.add(S[i-c+1:i+1])
        return lcs_set
    def topcs_seqrecords(sr1, sr2,x):
        naame=sr1.name+"__"+sr2.name
        if len(name) > 20:
            naame="__"+sr1[:7]
        lset=self.topcs( str(sr1.seq), str(sr2.seq),x)
        lsrecs=[]
        for ls in lset:
            lsrecs.append( SeqRecord( Seq(ls), id=naame, name="", description="" ) )
        return lsrecs
    def topcs_fragments(sr1, sr2,x):
        frags=[]
        tseqs=self.topcs_seqrecords(sr1,sr2, x)
        for i in tseqs:
            frags.append(str(tseqs.seq))
        return list(set(frags))
    def get_number_of_sequences(self):
        return len(self.farecords_)
    #x is the min fragment length
    def get_common_fragment_set_as_seqrecord_list(self, inda, indb, x=10):
        nid=self.farecords_[inda].id+"__"+self.farecords_[indb].id
        nname=self.farecords_[inda].name+"__"+self.farecords_[indb].name
        #return SeqRecord(seq=self.longestSubstring(self.farecords_[inda].seq, self.farecords_[indb].seq), id=nid, name=nname)
        seqsrno=1
        seqrecs=[]
        for e in self.topcs(str(self.farecords_[inda].seq), str(self.farecords_[indb].seq),x):
            sss=SeqRecord(Seq( e ,generic_dna),id=nid+"__"+str(seqsrno) , name=nname+"__"+str(seqsrno) )
            seqrecs.append(sss)
            seqsrno+=1
        return seqrecs
    def topcs_fragments(self, frag, seqrec, x=10):
        temprec=[]
        if type(seqrec) == "SeqRecord":
            temprec=self.topcs( frag, str(seqrec.seq), x )
        else:
            temprec=self.topcs(frag, seqrec, x)
        return list(set(temprec))
    def get_topcs_fragments(self, fraglist, seqreclist, x=10):
        tempfraglist=[]
        for s in seqreclist:
            for f in fraglist:
                tcf=self.topcs_fragments(f, s, x)
                for t in tcf:
                    if not t in tempfraglist:
                        tempfraglist.append(t)
        return tempfraglist
    def get_common_fragments(self,x):
        tempfrags=[]
        if self.get_number_of_sequences() == 1:
            tempfrags=[str(self.farecords_[0].seq)]
        elif self.get_number_of_sequences() > 1:
            tempfrags=self.topcs_fragments(self.farecords_[0], self.farecords_[1], x)
            tempfrags=self.get_topcs_fragments(tempfrags, self.farecords_[1:],x)
        return tempfrags
    def get_longest_common_fragment_set_as_seqrecord_list(self, inda, indb):
        nid=self.farecords_[inda].id+"__"+self.farecords_[indb].id
        nname=self.farecords_[inda].name+"__"+self.farecords_[indb].name
        #return SeqRecord(seq=self.longestSubstring(self.farecords_[inda].seq, self.farecords_[indb].seq), id=nid, name=nname)
        seqsrno=1
        seqrecs=[]
        for e in self.lcs(str(self.farecords_[inda].seq), str(self.farecords_[indb].seq)):
            sss=SeqRecord(Seq( e ,generic_dna),id=nid+"__"+str(seqsrno) , name=nname+"__"+str(seqsrno) )
            seqrecs.append(sss)
            seqsrno+=1
        return seqrecs
    #x is the min fragment length
    def write_common_framgments_as_fasta(self,x, outf):
        fo=open(outf, 'w')
        for i in range(1, self.get_number_of_sequences()):
            print("Computing for pair %d %d " %(0, i))
            ilist=self.get_common_fragment_set_as_seqrecord_list(0,i,x)
            print("LCSub generated for pair %d %d " %(0, i))
            print("Number of fragments found : %d" %(len(ilist)))
            for il in ilist:
                if len(il) > 0:
                    SeqIO.write(il, fo,"fasta")
        fo.close()
    def write_longest_common_framgments_as_fasta(self,outf):
        fo=open(outf, 'w')
        for i in range(1, self.get_number_of_sequences()):
            print("Computing for pair %d %d " %(0, i))
            ilist=self.get_longest_common_fragment_set_as_seqrecord_list(0,i)
            print("LCSub generated for pair %d %d " %(0, i))
            for il in ilist:
                if len(il) > 0:
                    SeqIO.write(il, fo,"fasta")
        fo.close()
def main():
    if len(sys.argv)<3:
        print("Usage: python3 LongestCommonFragment.py <inputfasta> <outputfasta>")
        sys.exit(0)
    else:
        l=LongestCommonFragment(sys.argv[1])
        #print("Object created!! Computing...")
        a=l.get_common_fragments(10)
        for i in a:
            print(a)
        #l.write_longest_common_framgments_as_fasta(sys.argv[2])
        #slist=l.get_longest_common_fragment_set_as_seqrecord_list(1,2)
        #for s in slist:
        #    print(s)

if __name__ == "__main__":
    main()

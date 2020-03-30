import sys

if len(sys.argv) < 2:
    print("usage:python3 agassv_to_tablerecords.py <aga_ssv_intermediate_file>")
    sys.exit(0)

class AGASSVRecord(object):
    def __init__(self, listrec):
        self.reference_=""
        self.query_=""
        self.genomealn_=[] #list of dictionaries //Every sequence may align multiple places in the genome or to multiple proteins or CDS
        self.cdsaln_=[] #list of dictionaries
        self.proteinaln_=[] #list of dictionaries
        if listrec[0].startswith("##reference"):
            self.reference_=listrec[0].split(":")[1]
            self.query_=listrec[1].split(":")[1]
            for li in range(2,len(listrec)):
                if listrec[li].startswith("#CDS"):
                    slrec=listrec[li][1:].strip().split(";")
                    cdsdict={}
                    for s in slrec:
                        if s != "":
                            ss=s.split(":")
                            cdsdict[ss[0]]=ss[1]
                    self.cdsaln_.append(cdsdict)
                elif listrec[li].startswith("#protein"):
                    slrec=listrec[li][1:].strip().split(";")
                    prodict={}
                    for s in slrec:
                        if s != "":
                            ss=s.split(":")
                            prodict[ss[0]]=ss[1]
                    self.proteinaln_.append(prodict)
                elif listrec[li].startswith("#genome"):
                    slrec=listrec[li][1:].strip().split(";")
                    gendict={}
                    for s in slrec:
                        if s != "":
                            ss=s.split(":")
                            gendict[ss[0]]=ss[1]
                    self.genomealn_.append(gendict)
        else:
            print("Format Error! Exiting code..")
            sys.exit(0)
        #Make it into a printable record string // right now only printing the first record, since that is relevant for current analysis //COV19
    def __str__(self):
        genomealndict=self.genomealn_[0]
        cdsalndict=self.cdsaln_[0]
        proteinalndict=self.proteinaln_[0]
        #print( len(self.proteinaln_))
        #print( len(self.cdsaln_))
        #print(len(self.genomealn_))
        if len(self.proteinaln_) == len(self.cdsaln_) and len(self.cdsaln_) == len(self.genomealn_):
            return str("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %( self.query_ , self.reference_ , genomealndict["genome"] ,  genomealndict["begin"], genomealndict["end"], genomealndict["matchpercent"], genomealndict["identitypercent"], cdsalndict["CDS"] ,   cdsalndict["begin"], cdsalndict["end"], cdsalndict["matchpercent"], cdsalndict["identitypercent"], proteinalndict["protein"], proteinalndict["begin"] , proteinalndict["end"], proteinalndict["matchpercent"], proteinalndict["identitypercent"]))
        else:
            llist=[ len(self.proteinaln_) , len(self.cdsaln_) , len(self.genomealn_)]
            maxv=max(llist)
            multirec=""
            for i in range(0, maxv):
                try:
                    genomealndict=self.genomealn_[i]
                except:
                    genomealndict=self.genomealn_[len(self.genomealn_)-1]
                try:
                    cdsalndict=self.cdsaln_[i]
                except:
                    cdsalndict=self.cdsaln_[len(self.cdsaln_)-1]
                try:
                    proteinalndict=self.proteinaln_[i]
                except:
                    proteinalndict=self.proteinaln_[len(self.proteinaln_)-1]
                multirec=multirec+str("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %( self.query_ , self.reference_ , genomealndict["genome"] ,  genomealndict["begin"], genomealndict["end"], genomealndict["matchpercent"], genomealndict["identitypercent"], cdsalndict["CDS"] ,   cdsalndict["begin"], cdsalndict["end"], cdsalndict["matchpercent"], cdsalndict["identitypercent"], proteinalndict["protein"], proteinalndict["begin"] , proteinalndict["end"], proteinalndict["matchpercent"], proteinalndict["identitypercent"]))+"\n"
            return multirec

def main():
    agassvrecs=[]
    fi=open(sys.argv[1])
    temprec=[]
    line1=fi.readline()
    if line1.startswith("##reference"):
        temprec.append(line1.strip())
    for i in fi:
        if i.startswith("##reference"):
            agassvrecs.append( AGASSVRecord(temprec) )
            del temprec[:]
            temprec.append(i.strip())
        else:
            temprec.append(i.strip())
    for a in agassvrecs:
        if a != None:
            print(str(a))
        else:
            pass
    fi.close()

if __name__=="__main__":
    main()

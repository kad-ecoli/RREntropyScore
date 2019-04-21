#!/usr/bin/env python
docstring='''
RR_ES_score.py T0859.seq.txt T0859.pdb T0859RR013_1_Lcont.json

    The three arguments are query sequence (fasta),
    native contact map (pdb or json), predicted contact map (rr or json)

Option:
    -infmt1={pdb,json} input format of native contact map
        pdb  - pdb structure
        json - CASP json format

    -infmt2={rr,gremlin,json} input format of predicted contact map
        rr      - CASP rr format
        gremlin - gremlin/ccmpred/deepcontact format
        json    - CASP json format

    -range={L,ML} sequence separation range
        L  - long range
        ML - medium and long range

    -atom={CA,CB} atom with which residue-residue distance is calculated
'''
import sys, os
import math
import json
import re
import numpy as np

def getSingleSeqFasta(input_fasta):
    ''' read single sequence fasta file and return the sequence '''
    fp=open(input_fasta,'rU')
    sequence=''.join([line.strip() for line in fp.read().splitlines(
        ) if not line.startswith('>')])
    fp.close()
    return sequence

class RRJsonManager:
    ''' 
    class to be rewritten in the future.
    It is recommended to read the original pdb file and RR map instead 
    '''
 
    # constructor
    def __init__(self):
        pass
  
    # read contacts in json format 
    def readConts(self, jsonfile):
        ''' 
        jsonfile - file in json format with data 
	containing contacts either in protein target or in prediction
	the data structures for these two cases are different
        '''
        fp=open(jsonfile)
        data = json.load(fp)
        fp.close()
        return data
  
    # reformat targets contacts eliminating the residues names
    def reformatTargetConts(self, data):
        result = {}
        for r1 in data:
            r1_ = r1
            lr2_ = []
            for r2 in data[r1]:
                r2_ = r2
                lr2_.append(r2_)
            result[r1_] = lr2_
        return result

    def readNativeConts(self, infile, sep_range, atom_sele="CB",
        dist_cutoff=8):
        ''' 
        infile      - input pdb file
        sep_range   - L for long range, ML for medium + long range
        atom_sele   - default is CA for gly, CB for all other atoms
        dist_cutoff - default is 8 angstrom
        '''
        result = dict()

        fp=sys.stdin
        if infile!='-':
            if infile.endswith(".gz"):
                fp=gzip.open(infile,'rU')
            else:
                fp=open(infile,'rU')
        struct=fp.read().split("ENDMDL")[0] # first model only
        fp.close()

        model=[r for r in struct.splitlines() if r.startswith("ATOM  ")]
        chain_id=[r[21] for r in model][0] # first chain
        chain=dict()
        for r in model:
            if r[21]!=chain_id:
                continue # first chain
            resName=r[17:20] 
            resSeq=r[22:26].strip()
            name=r[12:16].strip()
            x=float(r[30:38])
            y=float(r[38:46])
            z=float(r[46:54])
            if name==atom_sele or \
                (name=="CA" and not resSeq in chain): # CA if atom_sele is absent
                chain[resSeq]=(x,y,z)
        residues=sorted([k for k in chain]) # sorted list of residue index

        min_sep=sep_range
        if sep_range=="L":
            min_sep=24
        elif sep_range=="ML":
            min_sep=12
        for i in range(len(residues)-1):
            idx1=residues[i]
            x1,y1,z1=chain[idx1]
            if not idx1 in result:
                result[idx1]=[]
            for j in range(i+1,len(residues)):
                idx2=residues[j]
                if not idx2 in result:
                    result[idx2]=[]
                if abs(int(idx1)-int(idx2))<min_sep:
                    continue
                x2,y2,z2=chain[idx2]
                dx=x1-x2
                dy=y1-y2
                dz=z1-z2
                d=(dx*dx+dy*dy+dz*dz)**.5

                if d<=dist_cutoff:
                    result[idx1].append(idx2)
                    result[idx2].append(idx1)
        return result
  
    # filter prediction contacts by probability
    def filterRRcontacts(self, rr_data, targetConts, prob=0.0, rr_len = "FL"):
        if rr_len not in ("10", "L5", "L2", "FL"):
            rr_len = "FL"
        result = []
        for chain in rr_data:
            for cont in rr_data[chain][rr_len]:
                if float(cont[2]) >= prob:
                    try:
                        # check and adjust enumeration
                        if cont[1] in targetConts[cont[0]]:
                            result.append( (min(int(cont[0])-1, int(cont[1])-1), max(int(cont[0])-1, int(cont[1]))-1) )
                    except KeyError:
                        pass
            break
        return result

    def readRR(self, infile, prot_len, sep_range="L", offset=0, infmt="rr"):
        '''
        Read NN-BAYES or CASP RR format contact map. return them in a zipped 
        list with 3 fields for each residue pair. 1st field & 2nd filed are for 
        residue indices, and 3rd field is for euclidean distance.
        '''
        resi1=[] # residue index 1 list
        resi2=[] # residue index 2 list
        p=[] # cscore of contact prediction accuracy list
        fp=sys.stdin
        if infile!='-':
            if infile.endswith(".gz"):
                fp=gzip.open(infile,'rU')
            else:
                fp=open(infile,'rU')
        lines=fp.read().strip().splitlines()
        fp.close()
        pattern=re.compile('(^\d+\s+\d+\s+\d+\s+\d+\s+[-+.e\d]+)|(^\d+\s+\d+\s+[-+.e\d]+)|(^\d+\s+[A-Z]\s+\d+\s+[A-Z]\s+[-+.e\d]+\s+[-+.e\d]+)')

        min_sep=sep_range
        if sep_range=="L":
            min_sep=24
        elif sep_range=="ML":
            min_sep=12

        if infmt!="rr":
            for i,line in enumerate(lines):
                for j,cscore in enumerate(line.split()):
                    if abs(i-j)<min_sep:
                        continue
                    resi1.append(i+1+offset)
                    resi2.append(j+1+offset)
                    p.append(float(cscore))
        else:
            for line in lines:
                if not line.strip(): # skip empty lines
                    continue
                match_list=pattern.findall(line.strip())
                if not match_list:
                    continue
                line=[line for line in match_list[0] if line.strip()][0].split()
                if len(line)==6:
                    line=[line[0],line[2],line[5]]
                if not len(line) in (3,5):
                    continue

                resi_idx1=int(line[0])+offset # residue index 1
                resi_idx2=int(line[1])+offset # residue index 2
                cscore=float(line[-1]) # cscore for contact prediction
                if abs(resi_idx1-resi_idx2)<min_sep:
                    continue

                resi1.append(resi_idx1)
                resi2.append(resi_idx2)
                p.append(cscore)

        rr_data={'_':{'10':[], 'L5':[], 'L2':[], 'FL':[]}}

        cont_num=0
        for cscore,resi_idx1,resi_idx2 in sorted(zip(p,resi1,resi2),reverse=True):
            line=[str(resi_idx1),str(resi_idx2),str(cscore)]
            cont_num+=1
            if cont_num<=10:
                rr_data['_']['10'].append(line)
            if cont_num<=prot_len/2.:
                rr_data['_']['L2'].append(line)
                if cont_num<=prot_len/5.:
                    rr_data['_']['L5'].append(line)
            rr_data['_']['FL'].append(line)
        return rr_data

class RREntropyScore:
    # class constructor   
    def __init__(self):
        self.H_0 = {} # dictionary of uncertainty scores (key - length of a protein sequence) 

    # The method generates two matrices with upper and lower limits 
    # for a list of contacts
    def genUpperLowerLimitsMatrices(self, n, conts, l=3.2, u=8.0, L=3.2,
        U=None, neighbor_d = 3.8):
        '''
        n length of protein sequence, 
        conts - list of tuples (i,j) each representing the residue-residue contact
        l - lower limit for contact
        L - lower limit for non-contact
        u - upper limit for contact
        U - upper limit for non-contact 
        neighbor_d - distance between two neighbors; the default value corresponds to
        distance between two residues next to each other along the sequence
        '''
        if not U:
            # upper limit is set to the diameter of gyration calculated by Andriy's formula
            U = int(5.54*math.pow(n, 0.34))
        U_m = np.zeros((n,n))
        L_m = np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                if (i == j): # diagonal
                    U_m[i][j] = 0
                    L_m[i][j] = 0
                elif (abs(i - j) == 1): # neighbors
                    U_m[i][j] = neighbor_d
                    L_m[i][j] = neighbor_d
                else:
                    U_m[i][j] = U
                    L_m[i][j] = L

        for c in conts:
            i = c[0]
            j = c[1]
            L_m[i][j] = l
            U_m[i][j] = u
        return L_m, U_m

    def propogateTriangleInEqual(self, L_m, U_m):
        ''' 
        propagate the substitutions in upper and lower limit matrices
        A.W.M.Dress, T.F.Havel: Shortest-path problems and molecular conformation
	Disrt. Appl. Math. 19 (1988) 129-144
        '''
        n = len(L_m)
        for k in xrange(n):
            for i in xrange(n-1):
                for j in xrange(i+1, n):
                    if (U_m[i][j] > U_m[i][k]+U_m[k][j]):
                        U_m[i][j] = U_m[i][k]+U_m[k][j]
                        U_m[j][i] = U_m[i][j]
                    if (L_m[i][j] < L_m[i][k] - U_m[k][j]):
                        L_m[i][j] = L_m[i][k] - U_m[k][j]
                        L_m[j][i] = L_m[i][j] 
                    if (L_m[i][j] < L_m[j][k] - U_m[k][i]):
                        L_m[i][j] = L_m[j][k] - U_m[k][i]
                        L_m[j][i] = L_m[i][j]
                    if (L_m[i][j] > U_m[i][j]):
                        raise ValueError("Upper limit is less than lower limit for i=", i, ", j=", j)
        return L_m, U_m

    def calcUncertaintyScore(self, L_m, U_m):
        res = 0.0
        for i in xrange(len(L_m) - 1):
            for j in xrange(i + 1, len(L_m)):
                if (abs(i - j) <= 1 ): # ignore diagonal elements and "next door" neighbors
                    continue
                try:
                    res = res + math.log(U_m[i][j] - L_m[i][j])
                except ValueError:
                    print("Upper limit is less than lower limit for (", i,",",j,")")
                    return None
        return res    

    def calcScore(self, n, conts):
        '''
        Calculates the score based on Shannon's Information entropy:
        I = H|0 - H|r, structural information in r restraints (contacts)
          where H|0 - measure of uncertainty in structure without restraints
                H|r - measure of uncertainty in structure with restraints
        H = Sum(sum(log(U_ij - L_ij)))
        
        Ref: S. B. Nabuurs et al.: Quantitative Evaluation of Experimental
        NMR Restraints, J.Am.Chem.Soc., 2003, 125, 12026-12034        
        '''
        if not n in self.H_0.keys():
            # generate lower and upper limit matrices without restraints
            L_m0, U_m0 = self.genUpperLowerLimitsMatrices(n, [])
            L_m0, U_m0 = self.propogateTriangleInEqual(L_m0, U_m0)
            self.H_0[n] = self.calcUncertaintyScore(L_m0, U_m0)
        # generate lower and upper limit matrices with restraints
        L_m, U_m = self.genUpperLowerLimitsMatrices(n, conts)
        # propagate restraints through the matrices
        L_m, U_m = self.propogateTriangleInEqual(L_m, U_m)
        H_c = self.calcUncertaintyScore(L_m, U_m)
        score = 100.0*(self.H_0[n] - H_c)/(1.0*self.H_0[n])
        return score

if __name__=="__main__":
    infmt1="pdb"
    infmt2="rr"
    sep_range="L"
    atom_sele="CB"
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-infmt1="):
            infmt1=arg[len("-infmt1="):]
        elif arg.startswith("-infmt2="):
            infmt2=arg[len("-infmt2="):]
        elif arg.startswith("-range="):
            sep_range=arg[len("-range="):]
        elif arg.startswith("-atom="):
            atom_sele=arg[len("-atom="):]
        elif arg.startswith("-"):
            sys.stderr.write("ERROR! Unknown option %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)<3:
        sys.stderr.write(docstring)
        exit()
    if not infmt1 in {"pdb","json"}:
        sys.stderr.write("ERROR! Incorrect -infmt1=%s\n"%infmt1)
        exit()
    if not infmt2 in {"rr","gremlin","json"}:
        sys.stderr.write("ERROR! Incorrect -infmt1=%s\n"%infmt1)
        exit()

    input_fasta=argv[0]
    targetfile =argv[1]

    prot_len = len(getSingleSeqFasta(input_fasta))

    rrs = RREntropyScore() # invoke RREntropyScore class
    jm = RRJsonManager() # invoke RRJsonManager class

    if infmt1=="json":
        targetConts = jm.reformatTargetConts(jm.readConts(targetfile))
    else:
        targetConts = jm.readNativeConts(targetfile,sep_range,atom_sele)

    txt="Model,10p0,L5p0,L2p0,FLp0,10p05,L5p05,L2p05,FLp05\n"

    for modelfile in argv[2:]:
        #txt+='_'.join(modelfile.split('/')[-1].split('_')[:-1])+','
        #txt+=os.path.basename(modelfile).split('.')[0]+','
        txt+=modelfile+','

        if infmt2=="json":
            rr_data = jm.readConts(modelfile)
        else:
            rr_data = jm.readRR(modelfile, prot_len, sep_range)

        # loop over probabilities
        for _prob in [ 0.0, 0.5 ]:
            # loop over list sizes
            for _rr_len in [ "10", "L5", "L2", "FL" ]:
                rr_conts = jm.filterRRcontacts(rr_data, targetConts,
                    prob=_prob, rr_len=_rr_len)
                rrces=0.0
                if len(rr_conts):
                    rrces = rrs.calcScore(prot_len, rr_conts)
                txt+=str(rrces)
                if ( _rr_len == 'FL' and _prob == 0.5 ):
                    txt+="\n"
                else:
                    txt+=","

    print(txt.rstrip())

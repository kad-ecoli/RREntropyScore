#!/usr/bin/env python
docstring='''
RR_ES_score.py seqs/T0859.seq.txt                   \\
    input/T0859-D1_contList/T0859-D1_contactsL.json \\
    input/T0859-D1_contList/T0859RR013_1_Lcont.json

    The three arguments are query sequence (fasta),
    native contact map (json), predicted contact map (json)
'''
import sys
import math
import json
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
    class to be removed in the future
    it is recommended to read the original pdb file and RR map instead 
    '''
 
    # constructor
    def __init__(self):
        pass
  
    # read contacts in json format 
    def readConts(self, jsonfile):
        ''' 
        jsonfile - file in jason format with data 
	containing contacts either in protein target or in prediction
	the data structures for these two cases are different
        '''
        fp=open(jsonfile)
        data = json.load(fp)
        fp.close()
        return data
  
    # reformat targets contacts elemenibating the residues names
    def reformatTargetConts(self, data):
        result = {}
        for r1 in data:
            #r1_ = r1[3:] # obsolete format
            r1_ = r1
            lr2_ = []
            for r2 in data[r1]:
                #r2_ = r2[3:] # osolete format
                r2_ = r2
                lr2_.append(r2_)
            result[r1_] = lr2_
        return result
  
    # filter prediction contacts by probabilty
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
        propogate the substitutions in upper and lower limit matrices
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
        The method calculates the score based on Shannon's Information enthropy:
        I = H|0 - H|r, structural information in r restraints (contacts)
          where H|0 - measure of uncertainty in structure without restraints
                H|r - measure of uncertainty in structure with restraints
        H = Sum(sum(log(U_ij - L_ij)))
        
        Ref: S. B. Nabuurs et al.: Quantitative Evaluation of Experimenta; NMR Restraints, J.Am.Chem.Soc., 2003, 125, 12026-12034        
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
    if len(sys.argv)<3:
        sys.stderr.write(docstring)
        exit()

    input_fasta=sys.argv[1]
    targetfile =sys.argv[2]
    modelfile  =sys.argv[3]

    prot_len = len(getSingleSeqFasta(input_fasta))

    rrs = RREntropyScore() # invoke RREntropyScore class
    jm = RRJsonManager() # invoke RRJsonManager class

    targetConts = jm.reformatTargetConts(jm.readConts(targetfile))

    txt="Model,10p0,L5p0,L2p0,FLp0,10p05,L5p05,L2p05,FLp05\n"
    txt+='_'.join(modelfile.split('/')[-1].split('_')[:-1])+','

    rr_data = jm.readConts(modelfile)
    # loop over probabilities
    for _prob in [ 0.0, 0.5 ]:
        # loop over list sizes
        for _rr_len in [ "10", "L5", "L2", "FL" ]:
            rr_conts = jm.filterRRcontacts(rr_data, targetConts, prob=_prob, rr_len=_rr_len)
            rrces=0.0
            if len(rr_conts):
                rrces = rrs.calcScore(prot_len, rr_conts)
            txt+=str(rrces)
            if ( _rr_len == 'FL' and _prob == 0.5 ):
                txt+="\n"
            else:
                txt+=","
    print(txt.rstrip())

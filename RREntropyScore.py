
# import the required libraries
try: 
   import numpy as np
except ImportError:
   print "Failed to download librabry numpy"

try:
   import math
except ImportError:
   print "Failed to download librabry math"





class RREntropyScore:
  # class constructor   
  def __init__(self):
    self.H_0 = {} # dictionary of uncertainty scores (key - length of a protein sequence) 
    pass

  # The method generates two matrices with upper and lower limits 
  # for a list of contacts
  def genUpperLowerLimitsMatrices(self, n, conts, l=3.2, u=8.0, L=3.2, U=None, neighbor_d = 3.8):
    ''' n length of protein sequence, 
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
        else:
          if (abs(i - j) == 1): # neighbors
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
    ''' propogate the substitutions in upper and lower limit matrices
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
          print "Upper limit is less than lower limit for (", i,",",j,")"
          return None
    return res    

  def calcScore(self, n, conts):
    ''' The method calculates the score based on Shannon's Information enthropy:
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
   



# for testing/debugging
if __name__ == "__main__": 
  conts1 = [(10,50)]
  conts2 = [(10,50), (10,51), (11,50)]
  conts3 = [(10,25), (30,50), (70,90)]
  conts0 = []
  prot_len = 133
  rrs = RREntropyScore()
  print rrs.calcScore(prot_len, conts0)
  print rrs.calcScore(prot_len, conts1)
  print rrs.calcScore(prot_len, conts2)
  print rrs.calcScore(prot_len, conts3)


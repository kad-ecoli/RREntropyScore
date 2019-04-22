## RREntropyScore ##

The project is a re-implementation of entropy-based score (ES) for 
residue-residue contact predictions in CASP. Document for the original
version is available at [README](README).

#### ============================================== ####
The score calculates the relative drop of the entropy introduced 
by a set of distance constraints (residue-residue contacts) with 
the respect to the reference value of the entropy for the protein 
of a given length without constraints.
The score is calculated by formula:
```
ES = 100% (Entropy|0 - Entropy|C)/Entropy|0
```
where
```
Entropy|0 - the entropy value for the protein without constraints
Entropy|C - the entropy value given a set of constraints C
```
```
Entropy = AVERAGE over all pairs of residues (LOG(UpperLimit - LowerLimit))
```
where ``UpperLimit`` and ``LowerLimit`` stand for bounds for distances between
two residues. ``LowerLimit = 3.2``; ``UpperLimit = 8`` for contacts; 
``UpperLimit = DG - diameter of gyration`` for non-contacts. The diameter of 
gyration is estimated by formula ``DG=5.54L^0.34``
(``L - length of the protein sequence``).
(ref.: Nabuurs S.B et al: Quantitative Evaluation of Experimental NMR
Restraints, J. Am. Chem. Soc. 2003, 125, 12026-12034)

#### ============================================== ####
Instructions to run:
1. The program requires python module: numpy
2. Input: fasta format query sequence; pdb format native structure; CASP RR
   format contact map.
3. How to run from command line:
   ```bash
   python2 RR_ES_score.py seqs/T0859.seq.txt input/T0859-D1_contList/T0859.pdb input/T0859-D1_contList/T0859RR013_1.rr -range=L 
   ```
   Here, ``-range=L`` for long range contact; ``-range=ML`` for medium+long
   range contact. The above command will output
   ```
   Model,10p0,L5p0,L2p0,FLp0,10p05,L5p05,L2p05,FLp05
   input/T0859-D1_contList/T0859RR013_1.rr,0.0,0.109604354105,0.1566406265,1.27251629731,0.0,0.0,0.0,0.0
   ```
   Here, ``10``, ``L5``, ``F2``, and ``FL`` stand for top 10, top L/5,
   top L/2, and full contact map, respectively. ``p0`` and ``p05`` refer to
   ``cscore>=0`` and ``cscore>=0.5`` respectively.

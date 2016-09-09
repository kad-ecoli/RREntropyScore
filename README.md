# RREntropyScore

The project is related to calculation of entorpy-based score (ES) for 
residue-residue contact predictions in CASP.

==============================================
The score calculates the relative drop of the entropy introduced 
by a set of distance constraints (residue-residue contacts) with 
the respect to the reference value of the entropy for the protein 
of a given length without constraints.
The score is calculated by formula:
ES = 100% (Entropy|0 - Entropy|C)/Entropy|0
where Entropy|0 - the entropy value for the protein without constraints
Entropy|C - the entropy value given a set of constraints C
Entropy = AVERAGE over all pairs of residues (LOG(UpperLimit - LowerLimit)),
where UpperLimit and LowerLimit stand for bounds for distances between two residues, LowerLimit = 3.2Å
UpperLimit = 8Å for contacts; 
UpperLimit = DG - diameter of gyration for non-contacts,
The diameter of gyration is calculated by formula DG=5.54L^0.34 (L - length of the protein sequence).
(ref.: Nabuurs S.B et al: Quantitative Evaluation of Experimental NMR Restraints, J. Am. Chem. Soc. 2003, 125, 12026-12034)

==============================================


Instructions to run:
1. The program requires python module:
   1.1. numpy
2. Input: 
   2.1. lists of contacts in json format (Andriy's output)
        see examples in the directory ./input
   2.2. fasta sequence file for the CASP target in directory ./seqs
	see example
3. How to run from command line:
   python rrces.readconts.py <TARGET_NAME>
   e.g. python rrces.readconts.py T0859-D1
4. Output:
   The program generates files in directory ./output
   see example. 



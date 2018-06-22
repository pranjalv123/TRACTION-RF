# TRACTION-RF
Tree Refinement And CompleTION under Robinson-Foulds distances

TRACTION-RF takes as input a gene tree t (which may not be fully resolved) and a reference tree T (referred to below as the "backbone tree"), where the leaf set of t is a subset of the leaf set of T. It operates in two steps. In the first step, it refines t into a binary tree t' so as to minimize the total Robinson-Foulds (RF) distance to tree T (computed by restricting T to the leaf set of t). It then adds the missing species into tree t' to minimize the total RF distance to the tree T, using OCTAL (Christensen et al, Algorithms for Molecular Biology, 13:6, DOI: 10.1186/s13015-018-0124-5). The resultant two-step procedure provably produces a tree that minimizes the total RF distance to T. 

Current version: 1.0

Run TRACTION-RF as

    python traction.py -f -i <input gene trees> -b <backbone tree> -o <output> 

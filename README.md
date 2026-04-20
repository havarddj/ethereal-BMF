# Computations for "Ethereal Bianchi modular forms"

This directory contains computations relevant for the paper [CoD-J26]. 

- In [homology_torsion.m](./homology_torsion.m), there's code for computing the abelianization of $H_1$, which detects unliftable mod $p$ classes.
- In [builtin_utils.m](./builtin_utils.m) there are a couple of helper functions for computing and printing mod $p$ reductions of eigenvalues coming from the builtin Magma methods.
- In [write_forms.m](./write_forms.m) there's code for batch computing ethereal Bianchi classes.
- The directory [data](./data/) contains data computed using these methods.
- the directory [sage](./sage/) contains sage code for comparing eigenvalue systems with data from the LMFDB.

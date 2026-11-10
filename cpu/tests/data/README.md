# Small native-band regression fixture

`cdia_gamma.npz` contains 26 by 26 complex Hermitian H and matching nonorthogonal
S at Gamma, generated from this repository's `work/Cdia.dat` diamond example
(two carbon atoms, C5.0-s2p2d1 / C_CA19, LDA, 150 Ry, 7x7x7 SCF k grid).
It was exported with the unified CPU build based on develop `70021ccb` and
DFT_DATA19. `target` is the converged Fermi energy in Hartree.

The fixture reproduces an overly strict SLEPc BV inner-product roundoff check;
S is positive definite (eigenvalue range about 7.10e-4 to 8.58). The regression
compares exactly four requested eigenvalues against an independent full
**generalized** eigensolve and audits residuals and S orthogonality.
No user calculation or private path is embedded in the fixture.

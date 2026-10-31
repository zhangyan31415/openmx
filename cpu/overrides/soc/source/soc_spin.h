#ifndef SOC_SPIN_H
#define SOC_SPIN_H
/* Fields: UU.real, UU.imag, DD.real, DD.imag, UD.real, UD.imag,
 * DU.real, DU.imag. A directed atom/image block need not be Hermitian. */
int SocSpin_Transform(const double spin[8], double transform[64]);
int SocSpin_Reverse(int field, int *source, double *sign);
#endif

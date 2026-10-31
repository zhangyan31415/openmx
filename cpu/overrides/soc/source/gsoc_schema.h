#ifndef GSOC_SCHEMA_H
#define GSOC_SCHEMA_H
#include <math.h>
/* v1: proper R only. v2: physical R may be improper; metadata has already
 * constructed AO(R) times the SU(2) lift of the axial rotation det(R)*R.
 * The reader consumes the validated spinor CSC, never guesses a spin lift. */
static int gsoc_schema_version(int version){return version==1||version==2;}
static int gsoc_schema_determinant(int version,double determinant){
 if(!isfinite(determinant))return 0;
 if(version==1)return fabs(determinant-1.)<=1e-10;
 if(version==2)return fabs(fabs(determinant)-1.)<=1e-10;
 return 0;
}
#endif

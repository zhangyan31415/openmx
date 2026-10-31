#ifndef CPU_NSOC_IMAG_H
#define CPU_NSOC_IMAG_H
/* Pure policy: inactive imaginary history may alias current only in the
   separately validated owner real-Cluster fixed-geometry path. */
static inline int CpuNSOCImagSelect(int requested,int owner,int eligible,int history_slots)
{
  if (history_slots<1 || (requested!=0 && requested!=1)) return -1;
  if (!requested) return history_slots;
  if (history_slots<3) return -1;
  return owner && eligible ? 1 : -1;
}
int CpuNSOCImagSlots(void);
void CpuNSOCImagReport(unsigned long long local_payload_bytes);
#endif

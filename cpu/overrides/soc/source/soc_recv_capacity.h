#ifndef OPENMX_SOC_RECV_CAPACITY_H
#define OPENMX_SOC_RECV_CAPACITY_H
#include <stddef.h>
#include <stdint.h>
#include <limits.h>
/* Number of global indices in [0,end) owned by a block-cyclic coordinate.
 * Intermediates use size_t; no nb*np int multiplication is performed. */
static int SocRecvPrefix(int end,int nb,int coordinate,int np,size_t*out){
 if(end<0||nb<1||np<1||coordinate<0||coordinate>=np||!out)return -1;
 size_t full=(size_t)end/(size_t)nb,tail=(size_t)end%(size_t)nb;
 *out=(full/(size_t)np)*(size_t)nb+((size_t)coordinate<full%(size_t)np?(size_t)nb:0)+((size_t)coordinate==full%(size_t)np?tail:0);return 0;
}
/* Exact maximum receive count for the unchanged source producer: fields
 * values for each owned row and each owned column of a sector tile.
 * rank/row/col agreement is checked because legacy producer encodes rank as
 * row*np_cols+col. A different BLACS map was never valid for that producer. */
static int SocRecvMaximum(int n,int nb,int pr,int pc,int row,int col,int rank,
 int world,int tile,int fields,size_t*out){
 if(n<1||nb<1||pr<1||pc<1||(int64_t)pr*pc!=world||row<0||row>=pr||col<0||col>=pc||rank<0||rank>=world||(int64_t)row*pc+col!=rank||tile<1||fields<1||!out)return -1;
 size_t rows;if(SocRecvPrefix(n,nb,row,pr,&rows))return -1;size_t maximum=0;
 for(int begin=0;begin<n;){int width=n-begin<tile?n-begin:tile,end=begin+width;size_t lo,hi;
  if(SocRecvPrefix(begin,nb,col,pc,&lo)||SocRecvPrefix(end,nb,col,pc,&hi)||hi<lo)return -1;
  size_t cols=hi-lo;if(cols&&rows>SIZE_MAX/cols)return -1;size_t entries=rows*cols;
  if(entries>SIZE_MAX/(size_t)fields||entries*(size_t)fields>INT_MAX)return -1;
  size_t count=entries*(size_t)fields;if(count>maximum)maximum=count;begin=end;
 }
 if(maximum>SIZE_MAX/sizeof(double))return -1;*out=maximum;return 0;
}
static int SocRecvCountValid(size_t allocation,size_t required_max,int actual){
 return actual>=0&&(size_t)actual<=required_max&&(size_t)actual<=allocation;
}
static int SocRecvReadValid(size_t allocation,int actual,size_t offset,size_t width){
 if(actual<0||offset>allocation||width>allocation-offset)return 0;
 return offset<=(size_t)actual&&width<=(size_t)actual-offset;
}
#endif

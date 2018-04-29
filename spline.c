/*=======================================================================*/
#include <stdio.h>
#include <math.h>
#include "global_ex.h"
#include "global_var.h"
#include "ex_func.h"

int line1(pSPLINE pLine)
{
 float H[MAXNUM] = {0};
 float Fi[MAXNUM] = {0};
 float U[MAXNUM+1] = {0};
 float A[MAXNUM+1] = {0};
 float D[MAXNUM+1] = {0};
 float M[MAXNUM+1] = {0};
 float B[MAXNUM+1] = {0};
 float Y[MAXNUM+1] = {0};
 int i = 0;
 if((pLine->point_num < 3) || (pLine->point_num > MAXNUM + 1))
 {
  return -2;
 }
 for(i = 0;i <= pLine->point_num - 2;i++)
 {
  H[i] = pLine->x[i+1] - pLine->x[i];
  Fi[i] = (pLine->y[i+1] - pLine->y[i]) / H[i];
 }
for(i = 0;i <= pLine->point_num - 2;i++)
 {
  pLine->a3[i] = 0.0;
  pLine->a1[i] = Fi[i];
  pLine->b3[i] = 0.0;
  pLine->b1[i] = pLine->y[i];
 }
 return 1;

}


int Spline3(pSPLINE pLine)
{
 float H[MAXNUM] = {0};
 float Fi[MAXNUM] = {0};
 float U[MAXNUM+1] = {0};
 float A[MAXNUM+1] = {0};
 float D[MAXNUM+1] = {0};
 float M[MAXNUM+1] = {0};
 float B[MAXNUM+1] = {0};
 float Y[MAXNUM+1] = {0};
 int i = 0;
 if((pLine->point_num < 3) || (pLine->point_num > MAXNUM + 1))
 {
  return -2;
 }
 for(i = 0;i <= pLine->point_num - 2;i++)
 {
  H[i] = pLine->x[i+1] - pLine->x[i];
  Fi[i] = (pLine->y[i+1] - pLine->y[i]) / H[i];
 }
 for(i = 1;i <= pLine->point_num - 2;i++)
 {
  U[i] = H[i-1] / (H[i-1] + H[i]);
  A[i] = H[i] / (H[i-1] + H[i]);
  D[i] = 6 * (Fi[i] - Fi[i-1]) / (H[i-1] + H[i]);
 }
 
 U[i] = 0;
 A[0] = 0;
 pLine->begin_k2=0.0;
 pLine->end_k2=0.0;
 D[0] = 2 * pLine->begin_k2;
 D[i] = 2 * pLine->end_k2;

 B[0] = A[0] / 2;
 for(i = 1;i <= pLine->point_num - 2;i++)
 {
  B[i] = A[i] / (2 - U[i] * B[i-1]);
 }
 Y[0] = D[0] / 2;
 for(i = 1;i <= pLine->point_num - 1;i++)
 {
  Y[i] = (D[i] - U[i] * Y[i-1]) / (2 - U[i] * B[i-1]);
 }
 M[pLine->point_num - 1] = Y[pLine->point_num - 1];
 for(i = pLine->point_num - 1;i > 0;i--)
 {
  M[i-1] = Y[i-1] - B[i-1] * M[i];
 }

 for(i = 0;i <= pLine->point_num - 2;i++)
 {
  pLine->a3[i] = M[i] / (6 * H[i]);
  pLine->a1[i] = (pLine->y[i] - M[i] * H[i] * H[i] / 6) / H[i];
  pLine->b3[i] = M[i+1] / (6 * H[i]);
  pLine->b1[i] = (pLine->y[i+1] - M[i+1] * H[i] * H[i] / 6) /H[i];
 }
 return 1;
}

double func_spline3(double x,pSPLINE pLine)
{
    int i=MAXNUM+1;
    for(i=0;i<MAXNUM;i++){
        if( x >= pLine->x[i] && x<=pLine->x[i+1]){
        break;}
    }
        //printf("%d\n",i);
    return pLine->a3[i]*pow((pLine->x[i+1]-x),3)+pLine->a1[i]*(pLine->x[i+1]-x)+pLine->b3[i]*pow((x-pLine->x[i]),3)+pLine->b1[i]*(x-pLine->x[i]);


}

double func_line1(double x, pSPLINE pLine)
{
   int i=MAXNUM+1;
    for(i=0;i<MAXNUM;i++){
        if( x >= pLine->x[i] && x<=pLine->x[i+1]){
        break;}
    }

	return pLine->a1[i]*(x-pLine->x[i])+pLine->b1[i];
}


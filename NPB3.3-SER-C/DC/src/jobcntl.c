#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "adc.h"
#include "macrodef.h"

#ifdef UNIX
#include <fcntl.h>
#include <sys/file.h>
/*#include <sys/resource.h>*/
#include <unistd.h>
#endif

#include "adt_citerator.h"

#define USE_CITERATOR

uint32 NumberOfOnes(uint64 s);
void swap8(void *a);
void SetOneBit(uint64 *s, int32 pos){ uint64 ob = MLB; ob >>= pos; *s |= ob;}
void SetOneBit32(uint32 *s, uint32 pos){
   uint32 ob = 0x80000000;
   ob >>= pos;
   *s |= ob;
}
uint32 Mlo32(uint32 x){
   uint32 om = 0x80000000;
   uint32 i;
   uint32 k;

   for ( k = 0, i = 0; i < 32; i++ ) {
       if (om&x) break;
       om >>= 1;
       k++;
   }
   return(k);
}
int32 mro32(uint32 x){
   uint32 om = 0x00000001;
   uint32 i;
   uint32 k;

   for ( k = 32, i = 0; i < 32; i++ ) {
       if (om&x) break;
       om <<= 1;
       k--;
   }
   return(k);
}
uint32 setLeadingOnes32(uint32 n){
    int32 om = 0x80000000;
   uint32 x;
   uint32 i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
   FOR_START(x, cit1, 0, n, 1, cit_step_add, RND) {
   /*for ( x = 0, i = 0; i < n; i++ ) {*/
         x |= om;
         om >>= 1;
   }
   FOR_END(cit1);
#else
   for ( x = 0, i = 0; i < n; i++ ) {
         x |= om;
         om >>= 1;
   }
#endif // USE_CITERATOR
   return (x);
}
int32 DeleteOneFile(const char * file_name) {
#  ifdef WINNT
      return(remove(file_name));
#  else
      return(unlink(file_name));
#  endif
}
void WriteOne32Tuple(char * t, uint32 s, uint32 l, FILE * logf) {
  uint64 ob = MLB32;
  uint32 i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

  fprintf(logf, "\n %s", t);
#ifdef USE_CITERATOR
  FOR_START(i, cit1, 0, l, 1, cit_step_add, RND) {
  /*for ( i = 0; i < l; i++ ) {*/
    if (s&ob) fprintf(logf, "1"); else fprintf(logf, "0");
    ob >>= 1;
  }
  FOR_END(cit1);
#else
  for ( i = 0; i < l; i++ ) {
    if (s&ob) fprintf(logf, "1"); else fprintf(logf, "0");
    ob >>= 1;
  }
#endif // USE_CITERATOR
}
uint32 NumOfCombsFromNbyK( uint32 n, uint32 k ){
  uint32 l, combsNbyK;
  if ( k > n ) return 0;
  for(combsNbyK=1, l=1;l<=k;l++)combsNbyK = combsNbyK*(n-l+1)/l;
  return  combsNbyK;
}
void JobPoolUpdate(ADC_VIEW_CNTL *avp){
   uint32 l = avp->nv;
   uint32 k;

   k = avp->lpp[l].layerIndex + avp->lpp[l].layerCurrentPopulation;
   avp->jpp[k].grpb = avp->groupby;
   avp->jpp[k].nv = l;
   avp->jpp[k].nRows = avp->nViewRows;
   avp->jpp[k].viewOffset = avp->accViewFileOffset;
   avp->lpp[l].layerCurrentPopulation++;
}
int32 GetParent(ADC_VIEW_CNTL *avp, uint32 binRepTuple){
   uint32 level, levelPop, i;
   uint32 ig;
   uint32 igOfSmallestParent;
   uint32 igOfPrefixedParent;
   uint32 igOfSharedSortParent;
   uint32 spMinNumOfRows;
   uint32 pfMinNumOfRows;
   uint32 ssMinNumOfRows;
   uint32 tgrpb;
   uint32 pg;
   uint32 pfm;
   uint32 mlo = 0;
   uint32 lom;
   uint32 l = NumberOfOnes(binRepTuple);
   uint32 spFound;
   uint32 pfFound;
   uint32 ssFound;
   uint32 found;
   uint32 spFt;
   uint32 pfFt;
   uint32 ssFt;

   found = noneParent;
   pfm = setLeadingOnes32(mro32(avp->groupby));
   SetOneBit32(&mlo, Mlo32(avp->groupby));
   lom = setLeadingOnes32(Mlo32(avp->groupby));

   for(spFound=pfFound=ssFound=0, level=l;level<=avp->nTopDims;level++){
      levelPop = avp->lpp[level].layerCurrentPopulation;

      if(levelPop != 0);
      {
           for ( spFt = pfFt = ssFt = 1, ig = avp->lpp[level].layerIndex,
                 i = 0; i < levelPop; i++ )
           {
               tgrpb = avp->jpp[ig].grpb;
               if ( (avp->groupby & tgrpb) == avp->groupby ) {
                  spFound = 1;
                  if (spFt) { spMinNumOfRows = avp->jpp[ig].nRows;
                              igOfSmallestParent = ig; spFt = 0; }
                  else   if ( spMinNumOfRows > avp->jpp[ig].nRows )
                            { spMinNumOfRows = avp->jpp[ig].nRows;
                              igOfSmallestParent = ig; }

				  pg = tgrpb & pfm;
				  if (pg == binRepTuple) {
                     pfFound = 1;
                     if (pfFt) { pfMinNumOfRows = avp->jpp[ig].nRows;
                                 igOfPrefixedParent = ig; pfFt = 0; }
                     else   if ( pfMinNumOfRows > avp->jpp[ig].nRows)
                               { pfMinNumOfRows = avp->jpp[ig].nRows;
                                 igOfPrefixedParent = ig; }
				  }

				  if ( (tgrpb & mlo) && !(tgrpb & lom)) {
                     ssFound = 1;
                     if (ssFt) { ssMinNumOfRows = avp->jpp[ig].nRows;
                                 igOfSharedSortParent = ig; ssFt = 0; }
                     else   if ( ssMinNumOfRows > avp->jpp[ig].nRows)
                               { ssMinNumOfRows = avp->jpp[ig].nRows;
                                 igOfSharedSortParent = ig; }
				  }
               }
               ig++;
           }
      }
      if (pfFound) found = prefixedParent;
      else if (ssFound) found = sharedSortParent;
           else if (spFound) found = smallestParent;

      switch(found){
         case prefixedParent:
           avp->smallestParentLevel = level;
           avp->viewOffset      = avp->jpp[igOfPrefixedParent].viewOffset;
           avp->nParentViewRows = avp->jpp[igOfPrefixedParent].nRows;
           avp->parBinRepTuple  = avp->jpp[igOfPrefixedParent].grpb;
           break;
         case sharedSortParent:
           avp->smallestParentLevel = level;
           avp->viewOffset	    = avp->jpp[igOfSharedSortParent].viewOffset;
           avp->nParentViewRows = avp->jpp[igOfSharedSortParent].nRows;
           avp->parBinRepTuple  = avp->jpp[igOfSharedSortParent].grpb;
           break;
         case smallestParent:
           avp->smallestParentLevel = level;
           avp->viewOffset	    = avp->jpp[igOfSmallestParent].viewOffset;
           avp->nParentViewRows = avp->jpp[igOfSmallestParent].nRows;
           avp->parBinRepTuple  = avp->jpp[igOfSmallestParent].grpb;
           break;
         default: break;
      }
      if(   found == prefixedParent
         || found == sharedSortParent
	 || found == smallestParent) break;
   }
  return found;
}
uint32 GetSmallestParent(ADC_VIEW_CNTL *avp, uint32 binRepTuple){
   uint32 found, level, levelPop, i, ig, igOfSmallestParent;
   uint32 minNumOfRows;
   uint32 tgrpb;
   uint32 ft;
   uint32 l = NumberOfOnes(binRepTuple);

   for(found=0, level=l; level<=avp->nTopDims;level++){
      levelPop = avp->lpp[level].layerCurrentPopulation;
      if(levelPop){
        for(ft=1, ig=avp->lpp[level].layerIndex, i=0;i<levelPop;i++){
          tgrpb = avp->jpp[ig].grpb;
          if ( (avp->groupby & tgrpb) == avp->groupby ) {
            found = 1;
            if(ft){
	      minNumOfRows=avp->jpp[ig].nRows;
	      igOfSmallestParent = ig;
	      ft = 0;
	    }else if(minNumOfRows > avp->jpp[ig].nRows){
	      minNumOfRows = avp->jpp[ig].nRows;
	      igOfSmallestParent = ig;
	    }
          }
          ig++;
        }
      }
      if( found ){
         avp->smallestParentLevel = level;
         avp->viewOffset = avp->jpp[igOfSmallestParent].viewOffset;
         avp->nParentViewRows = avp->jpp[igOfSmallestParent].nRows;
         avp->parBinRepTuple = avp->jpp[igOfSmallestParent].grpb;
         break;
      }
   }
   return found;
}
int32 GetPrefixedParent(ADC_VIEW_CNTL *avp, uint32 binRepTuple){
   uint32 found, level, levelPop, i, ig, igOfSmallestParent;
   uint32 minNumOfRows;
   uint32 tgrpb;
   uint32 ft;
   uint32 pg, tm;
   uint32 l = NumberOfOnes(binRepTuple);

   tm = setLeadingOnes32(mro32(avp->groupby));

   for(found=0, level=l; level<=avp->nTopDims; level++){
      levelPop = avp->lpp[level].layerCurrentPopulation;

      if (levelPop != 0);
      {
           for(ft = 1, ig = avp->lpp[level].layerIndex,
                i = 0; i < levelPop; i++ ) {
               tgrpb = avp->jpp[ig].grpb;
               if ( (avp->groupby & tgrpb) == avp->groupby ) {
				  pg = tgrpb & tm;
				  if (pg == binRepTuple) {
                     found = 1;
                     if (ft) { minNumOfRows = avp->jpp[ig].nRows;
                               igOfSmallestParent = ig; ft = 0; }
                     else if ( minNumOfRows > avp->jpp[ig].nRows)
                             { minNumOfRows = avp->jpp[ig].nRows;
                               igOfSmallestParent = ig; }
				  }
               }
               ig++;
           }
      }
      if ( found ) {
         avp->smallestParentLevel = level;
         avp->viewOffset = avp->jpp[igOfSmallestParent].viewOffset;
         avp->nParentViewRows = avp->jpp[igOfSmallestParent].nRows;
         avp->parBinRepTuple = avp->jpp[igOfSmallestParent].grpb;
         break;
      }
   }
  return found;
}
void JobPoolInit(JOB_POOL *jpp, uint32 n, uint32 nd){
  uint32 i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, 0, n, 1, cit_step_add, RND) {
  /*for ( i = 0; i < n; i++ ) {*/
      jpp[i].grpb = 0;
	  jpp[i].nv = 0;
      jpp[i].nRows = 0;
      jpp[i].viewOffset = 0;
  }
  FOR_END(cit1);
#else
  for ( i = 0; i < n; i++ ) {
      jpp[i].grpb = 0;
	  jpp[i].nv = 0;
      jpp[i].nRows = 0;
      jpp[i].viewOffset = 0;
  }
#endif // USE_CITERATOR
}
void WriteOne64Tuple(char * t, uint64 s, uint32 l, FILE * logf){
   uint64 ob = MLB;
   uint32 i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

   fprintf(logf, "\n %s", t);
#ifdef USE_CITERATOR
   FOR_START(i, cit1, 0, l, 1, cit_step_add, RND) {
   /*for ( i = 0; i < l; i++ ) {*/
      if (s&ob) fprintf(logf, "1"); else fprintf(logf, "0");
      ob >>= 1;
   }
   FOR_END(cit1);
#else
   for ( i = 0; i < l; i++ ) {
      if (s&ob) fprintf(logf, "1"); else fprintf(logf, "0");
      ob >>= 1;
   }
#endif // USE_CITERATOR
}
uint32 NumberOfOnes(uint64 s){
   uint64 ob = MLB;
   uint32 i;
   uint32 nOnes;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
   FOR_START(nOnes, cit1, 0, 64, 1, cit_step_add, RND) {
   /*for ( nOnes = 0, i = 0; i < 64; i++ ) {*/
      if (s&ob) nOnes++;
      ob >>= 1;
   }
   FOR_END(cit1);
#else
   for ( nOnes = 0, i = 0; i < 64; i++ ) {
      if (s&ob) nOnes++;
      ob >>= 1;
   }
#endif // USE_CITERATOR
   return nOnes;
}
void GetRegTupleFromBin64(
           uint64 binRepTuple,
	       uint32 *selTuple,
	       uint32 numDims,
	       uint32 *numOfUnits){
   uint64 oc = MLB;
   uint32 i;
   uint32 j;
#ifdef USE_CITERATOR
   struct cit_data *cit1;
#endif // USE_CITERATOR

   *numOfUnits = 0;
#ifdef USE_CITERATOR
   j = 0;
   FOR_START(i, cit1, 0, numDims, 1, cit_step_add, FWD) {
   /*for( j = 0, i = 0; i < numDims; i++ ) {*/
     if (binRepTuple & oc) { selTuple[j++] = i+1; (*numOfUnits)++;}
     oc >>= 1;
   }
   FOR_END(cit1);
#else
   for( j = 0, i = 0; i < numDims; i++ ) {
     if (binRepTuple & oc) { selTuple[j++] = i+1; (*numOfUnits)++;}
     oc >>= 1;
   }
#endif // USE_CITERATOR
}
void getRegTupleFromBin32(
           uint32 binRepTuple,
	       uint32 *selTuple,
	       uint32 numDims,
	       uint32 *numOfUnits){
   uint32 oc = MLB32;
   uint32 i;
   uint32 j;
#ifdef USE_CITERATOR
   struct cit_data *cit1;
#endif // USE_CITERATOR

   *numOfUnits = 0;
#ifdef USE_CITERATOR
   j = 0;
   FOR_START(i, cit1, 0, numDims, 1, cit_step_add, RND) {
   /*for( j = 0, i = 0; i < numDims; i++ ) {*/
     if (binRepTuple & oc) { selTuple[j++] = i+1; (*numOfUnits)++;}
     oc >>= 1;
   }
   FOR_END(cit1);
#else
   for( j = 0, i = 0; i < numDims; i++ ) {
     if (binRepTuple & oc) { selTuple[j++] = i+1; (*numOfUnits)++;}
     oc >>= 1;
   }
#endif // USE_CITERATOR
}
void GetRegTupleFromParent(
               uint64 bin64RepTuple,
               uint32 bin32RepTuple,
	       uint32 *selTuple,
	       uint32 nd){
   uint32 oc = MLB32;
   uint32 i, j, k;
   uint32 ut32;
#ifdef USE_CITERATOR
   struct cit_data *cit1;
#endif // USE_CITERATOR

   ut32 = (uint32)(bin64RepTuple>>(64-nd));
   ut32 <<= (32-nd);

#ifdef USE_CITERATOR
   j = k = 0;
   FOR_START(i, cit1, 0, nd, 1, cit_step_add, RND) {
   /*for ( j = 0, k = 0, i = 0; i < nd; i++ ) {*/
     if (bin32RepTuple & oc) k++;
     if (bin32RepTuple & oc && ut32 & oc) selTuple[j++] = k;
     oc >>= 1;
   }
   FOR_END(cit1);
#else
   for ( j = 0, k = 0, i = 0; i < nd; i++ ) {
     if (bin32RepTuple & oc) k++;
     if (bin32RepTuple & oc && ut32 & oc) selTuple[j++] = k;
     oc >>= 1;
   }
#endif // USE_CITERATOR
}
void CreateBinTuple(uint64 *binRepTuple, uint32 *selTuple, uint32 numDims){
   uint32 i;
#ifdef USE_CITERATOR
   struct cit_data *cit1;
#endif // USE_CITERATOR

   *binRepTuple = 0;
#ifdef USE_CITERATOR
   FOR_START(i, cit1, 0, numDims, 1, cit_step_add, RND) {
   /*for(i = 0; i < numDims; i++ ){*/
     SetOneBit( binRepTuple, selTuple[i]-1 );
   }
   FOR_END(cit1);
#else
   for(i = 0; i < numDims; i++ ){
     SetOneBit( binRepTuple, selTuple[i]-1 );
   }
#endif // USE_CITERATOR
}
void d32v( char * t, uint32 *v, uint32 n){
   uint32 i;

   fprintf(stderr,"\n%s ", t);
   for ( i = 0; i < n; i++ ) fprintf(stderr," %d", v[i]);
}
void WriteOne64Tuple(char * t, uint64 s, uint32 l, FILE * logf);
int32 Comp8gbuf(const void *a, const void *b){
   if ( a < b ) return -1;
   else if (a > b) return 1;
   else return 0;
}
void restore(TUPLE_VIEWSIZE x[], uint32 f, uint32 l ){
   uint32 j, m, tj, mm1, jm1, hl;
   uint64 iW;
   uint64 iW64;

   j = f;
   hl = l>>1;
   while( j <= hl ) {
      tj = j*2;
      if (tj < l && x[tj-1].viewsize < x[tj].viewsize) m = tj+1;
      else m = tj;
      mm1 = m - 1;
      jm1 = j - 1;
      if ( x[mm1].viewsize > x[jm1].viewsize ) {
         iW = x[mm1].viewsize;
	 x[mm1].viewsize = x[jm1].viewsize;
	 x[jm1].viewsize = iW;
         iW64 = x[mm1].tuple;
	 x[mm1].tuple = x[jm1].tuple;
	 x[jm1].tuple = iW64;
         j = m;
      }else j = l;
   }
}
void vszsort( TUPLE_VIEWSIZE x[], uint32 n){
  int32 i, im1;
  uint64 iW;
  uint64 iW64;
#ifdef USE_CITERATOR
    struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, n>>1, 0, -1, cit_step_add, RND) {
  /*for ( i = n>>1; i >= 1; i-- ) {*/
    restore( x, i, n );
  }
  FOR_END(cit1);
#else
  for ( i = n>>1; i >= 1; i-- ) restore( x, i, n );
#endif // USE_CITERATOR
#ifdef USE_CITERATOR
  FOR_START(i, cit1, n, 1, -1, cit_step_add, RND) {
  /*for ( i = n; i >= 2; i-- ) {*/
     im1 = i - 1;
     iW = x[0].viewsize; x[0].viewsize = x[im1].viewsize; x[im1].viewsize = iW;
     iW64 = x[0].tuple; x[0].tuple = x[im1].tuple; x[im1].tuple = iW64;
     restore( x, 1, im1);
  }
  FOR_END(cit1);
#else
  for ( i = n; i >= 2; i-- ) {
     im1 = i - 1;
     iW = x[0].viewsize; x[0].viewsize = x[im1].viewsize; x[im1].viewsize = iW;
     iW64 = x[0].tuple; x[0].tuple = x[im1].tuple; x[im1].tuple = iW64;
     restore( x, 1, im1);
  }
#endif // USE_CITERATOR
}
uint32 countTupleOnes(uint64 binRepTuple, uint32 numDims){
  uint32 i, cnt = 0;
  uint64 ob = 0x0000000000000001;

  for(i = 0; i < numDims; i++ ){
    if ( binRepTuple&ob) cnt++;
    ob <<= 1;
  }
  return cnt;
}
void restoreo( TUPLE_ONES x[], uint32 f, uint32 l ){
   uint32 j, m, tj, mm1, jm1, hl;
   uint32 iW;
   uint64 iW64;

   j = f;
   hl = l>>1;
   while( j <= hl ) {
      tj = j*2;
      if (tj < l && x[tj-1].nOnes < x[tj].nOnes) m = tj+1;
      else m = tj;
      mm1 = m - 1; jm1 = j - 1;
      if ( x[mm1].nOnes > x[jm1].nOnes ){
         iW = x[mm1].nOnes;
	     x[mm1].nOnes = x[jm1].nOnes;
	     x[jm1].nOnes = iW;
         iW64 = x[mm1].tuple;
	     x[mm1].tuple = x[jm1].tuple;
	     x[jm1].tuple = iW64;
         j = m;
      }else j = l;
   }
}
void onessort( TUPLE_ONES x[], uint32 n){
   int32 i, im1;
  uint32 iW;
  uint64 iW64;
#ifdef USE_CITERATOR
    struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, n>>1, 0, -1, cit_step_add, RND) {
  /*for ( i = n>>1; i >= 1; i-- ) {*/
    restoreo( x, i, n );
  }
  FOR_END(cit1);
#else
  for ( i = n>>1; i >= 1; i-- ) restoreo( x, i, n );
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, n, 1, -1, cit_step_add, RND) {
  /*for ( i = n; i >= 2; i-- ) {*/
     im1 = i - 1;
     iW = x[0].nOnes;
     x[0].nOnes = x[im1].nOnes;
     x[im1].nOnes = iW;
     iW64 = x[0].tuple;
     x[0].tuple = x[im1].tuple;
     x[im1].tuple = iW64;
     restoreo( x, 1, im1);
  }
  FOR_END(cit1);
#else
  for ( i = n; i >= 2; i-- ) {
     im1 = i - 1;
     iW = x[0].nOnes;
     x[0].nOnes = x[im1].nOnes;
     x[im1].nOnes = iW;
     iW64 = x[0].tuple;
     x[0].tuple = x[im1].tuple;
     x[im1].tuple = iW64;
     restoreo( x, 1, im1);
  }
#endif // USE_CITERATOR
}
uint32 MultiFileProcJobs( TUPLE_VIEWSIZE *tuplesAndSizes,
		                          uint32 nViews,
                           ADC_VIEW_CNTL *avp ){
   uint32 i;
    int32 ii; /* it should be int */
   uint32 j;
   uint32 pn;
   uint32 direction = 0;
   uint32 dChange = 0;
   uint32 gbi;
   uint32 maxn;
   uint64 *gbuf;
   uint64      vszs[MAX_NUMBER_OF_TASKS];
   uint32 nGroupbys[MAX_NUMBER_OF_TASKS];
   TUPLE_ONES *toptr;
#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2;
#endif // USE_CITERATOR

   gbuf = (uint64*) &avp->memPool[0];

#ifdef USE_CITERATOR
   FOR_START(i, cit1, 0, avp->nTasks, 1, cit_step_add, RND) {
   /*for(i = 0; i < avp->nTasks; i++ ){*/
     nGroupbys[i] = 0; vszs[i] = 0;
   }
#else
   for(i = 0; i < avp->nTasks; i++ ){ nGroupbys[i] = 0; vszs[i] = 0; }
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
   pn=0;
   gbi=0;
   FOR_START(ii, cit1, nViews-1, -1, -1, cit_step_add, RND) {
   /*for(pn = 0, gbi = 0, ii = nViews-1; ii >= 0; ii-- ){*/
     if(pn == avp->taskNumber) gbuf[gbi++]=tuplesAndSizes[ii].tuple;
     nGroupbys[pn]++;
     vszs[pn] += tuplesAndSizes[ii].viewsize;
     if(direction == 0 && pn == avp->nTasks-1 ) {
       direction = 1;
       dChange = 1;
     }
     if(direction == 1 && pn == 0 ){
       direction = 0;
       dChange = 1;
     }
     if (!dChange){ if (direction) pn--; else pn++;}
     dChange = 0;
   }
   FOR_END(cit1);
#else
   for(pn = 0, gbi = 0, ii = nViews-1; ii >= 0; ii-- ){
     if(pn == avp->taskNumber) gbuf[gbi++]=tuplesAndSizes[ii].tuple;
     nGroupbys[pn]++;
     vszs[pn] += tuplesAndSizes[ii].viewsize;
     if(direction == 0 && pn == avp->nTasks-1 ) {
       direction = 1;
       dChange = 1;
     }
     if(direction == 1 && pn == 0 ){
       direction = 0;
       dChange = 1;
     }
     if (!dChange){ if (direction) pn--; else pn++;}
     dChange = 0;
   }
#endif // USE_CITERATOR
   for(maxn = 0, i = 0; i < avp->nTasks; i++)
     if (nGroupbys[i] > maxn) maxn = nGroupbys[i];

   toptr = (TUPLE_ONES*) malloc(sizeof(TUPLE_ONES)*maxn);
   if(!toptr) return 1;

#ifdef USE_CITERATOR
   FOR_START(i, cit1, 0, avp->nTasks, 1, cit_step_add, RND) {
   /*for(i = 0; i < avp->nTasks; i++ ){*/
     if(i == avp->taskNumber){
       FOR_START(j, cit2, 0, nGroupbys[i], 1, cit_step_add, RND) {
       /*for(j = 0; j < nGroupbys[i]; j++ ){*/
         toptr[j].tuple = gbuf[j];
         toptr[j].nOnes  = countTupleOnes(gbuf[j], avp->nTopDims);
       }
       FOR_END(cit2);
       qsort((void*)gbuf,  nGroupbys[i], 8, Comp8gbuf );
       onessort(toptr, nGroupbys[i]);

       FOR_START(j, cit2, 0, nGroupbys[i], 1, cit_step_add, RND) {
       /*for(j = 0; j < nGroupbys[i]; j++){*/
         toptr[nGroupbys[i]-1-j].tuple <<= (64-avp->nTopDims);
         swap8(&toptr[nGroupbys[i]-1-j].tuple);
         fwrite(&toptr[nGroupbys[i]-1-j].tuple, 8, 1, avp->groupbyFile);
       }
       FOR_END(cit2);
     }
   }
   FOR_END(cit1);
#else
   for(i = 0; i < avp->nTasks; i++ ){
     if(i == avp->taskNumber){
       for(j = 0; j < nGroupbys[i]; j++ ){
         toptr[j].tuple = gbuf[j];
         toptr[j].nOnes  = countTupleOnes(gbuf[j], avp->nTopDims);
       }
       qsort((void*)gbuf,  nGroupbys[i], 8, Comp8gbuf );
       onessort(toptr, nGroupbys[i]);

       for(j = 0; j < nGroupbys[i]; j++){
         toptr[nGroupbys[i]-1-j].tuple <<= (64-avp->nTopDims);
         swap8(&toptr[nGroupbys[i]-1-j].tuple);
         fwrite(&toptr[nGroupbys[i]-1-j].tuple, 8, 1, avp->groupbyFile);
       }
     }
   }
#endif // USE_CITERATOR
   FSEEK(avp->groupbyFile, 0L, SEEK_SET);
   if (toptr) free(toptr);
   return 0;
}
int32 PartitionCube(ADC_VIEW_CNTL *avp){
    TUPLE_VIEWSIZE *tuplesAndSizes;
    uint32 it = 0;
    uint64 sz;
    uint32 sel[64];
    uint32 k;
    uint64 tx;
    uint32 i;
      char inps[256];
#ifdef USE_CITERATOR
    struct cit_data *cit1;
#endif // USE_CITERATOR

    tuplesAndSizes =
       (TUPLE_VIEWSIZE*) malloc(avp->nViewLimit*sizeof(TUPLE_VIEWSIZE));
    if(tuplesAndSizes == NULL){
       fprintf(stderr," PartitionCube(): memory allocation failure'\n");
       return ADC_MEMORY_ALLOCATION_FAILURE;
    }
    k = 0;
    while( fscanf(avp->adcViewSizesFile, "%s", inps) != EOF ){
       if( strcmp(inps, "Selection:") == 0 ) {
         while ( fscanf(avp->adcViewSizesFile, "%s", inps)) {
           if ( strcmp(inps, "View") == 0 ) break;
           sel[k++] = atoi(inps);
         }
       }
       if( strcmp(inps, "Size:") == 0 ){
         fscanf(avp->adcViewSizesFile, "%s", inps);
         sz = atoi(inps);
         CreateBinTuple(&tx, sel, k);
         if (sz > avp->nInputRecs) sz = avp->nInputRecs;
         tuplesAndSizes[it].viewsize = sz;
         tuplesAndSizes[it].tuple = tx;
         it++;
         k = 0;
       }
    }
    vszsort(tuplesAndSizes, it);
#ifdef USE_CITERATOR
    FOR_START(i, cit1, 0, it, 1, cit_step_add, RND) {
    /*for( i = 0; i < it; i++){*/
        tuplesAndSizes[i].tuple >>= (64-avp->nTopDims);
    }
    FOR_END(cit1);
#else
    for( i = 0; i < it; i++){
        tuplesAndSizes[i].tuple >>= (64-avp->nTopDims);
    }
#endif // USE_CITERATOR
    if(MultiFileProcJobs( tuplesAndSizes, it, avp )){
       fprintf(stderr, "MultiFileProcJobs() is failed \n");
       fprintf(avp->logf, "MultiFileProcJobs() is failed.\n");
       fflush(avp->logf);
       return 1;
    }
    FSEEK(avp->adcViewSizesFile, 0L, SEEK_SET);
    free(tuplesAndSizes);
    return 0;
}

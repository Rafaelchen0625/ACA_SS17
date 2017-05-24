#include <immintrin.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#define IDCT_SIZE         16
#define ITERATIONS        1000000
#define MAX_NEG_CROP      1024

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

static const short g_aiT16[16][16] =
{
  { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64},
  { 90, 87, 80, 70, 57, 43, 25,  9, -9,-25,-43,-57,-70,-80,-87,-90},
  { 89, 75, 50, 18,-18,-50,-75,-89,-89,-75,-50,-18, 18, 50, 75, 89},
  { 87, 57,  9,-43,-80,-90,-70,-25, 25, 70, 90, 80, 43, -9,-57,-87},
  { 83, 36,-36,-83,-83,-36, 36, 83, 83, 36,-36,-83,-83,-36, 36, 83},
  { 80,  9,-70,-87,-25, 57, 90, 43,-43,-90,-57, 25, 87, 70, -9,-80},
  { 75,-18,-89,-50, 50, 89, 18,-75,-75, 18, 89, 50,-50,-89,-18, 75},
  { 70,-43,-87,  9, 90, 25,-80,-57, 57, 80,-25,-90, -9, 87, 43,-70},
  { 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64, 64,-64,-64, 64},
  { 57,-80,-25, 90, -9,-87, 43, 70,-70,-43, 87,  9,-90, 25, 80,-57},
  { 50,-89, 18, 75,-75,-18, 89,-50,-50, 89,-18,-75, 75, 18,-89, 50},
  { 43,-90, 57, 25,-87, 70,  9,-80, 80, -9,-70, 87,-25,-57, 90,-43},
  { 36,-83, 83,-36,-36, 83,-83, 36, 36,-83, 83,-36,-36, 83,-83, 36},
  { 25,-70, 90,-80, 43,  9,-57, 87,-87, 57, -9,-43, 80,-90, 70,-25},
  { 18,-50, 75,-89, 89,-75, 50,-18,-18, 50,-75, 89,-89, 75,-50, 18},
  {  9,-25, 43,-57, 70,-80, 87,-90, 90,-87, 80,-70, 57,-43, 25, -9}
};

static int64_t diff(struct timespec start, struct timespec end)
{
    struct timespec temp;
    int64_t d;
    if ((end.tv_nsec-start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    d = temp.tv_sec*1000000000+temp.tv_nsec;
    return d;
}

static void compare_results(short *ref, short *res, const char *msg)
{
    int correct =1;

    printf("Comparing %s\n",msg);
    for(int j=0; j<IDCT_SIZE; j++)  {
        for(int i=0; i<IDCT_SIZE; i++){
            if(ref[j*IDCT_SIZE+i] != res[j*IDCT_SIZE+i]){
                correct=0;
                printf("failed at %d,%d\t ref=%d, res=%d\n ", i, j, ref[j*IDCT_SIZE+i],res[j*IDCT_SIZE+i]);
            }
        }
    }
    if (correct){
        printf("correct\n\n");
    }
}

// this function is for timing, do not change anything here
static void benchmark( void (*idct16)(short *, short *), short *input, short *output, const char *version )
{
    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME,&start);

    for(int i=0;i<ITERATIONS;i++)
        idct16(input, output);

    clock_gettime(CLOCK_REALTIME,&end);
    double avg = (double) diff(start,end)/ITERATIONS;
    printf("%10s:\t %.3f ns\n", version, avg);
}

//scalar code for the inverse transform
static void partialButterflyInverse16(short *src, short *dst, int shift)
{
  int E[8],O[8];
  int EE[4],EO[4];
  int EEE[2],EEO[2];
  int add = 1<<(shift-1);

  //new Matrix: each row g = g_aiT16[i]*src[i*16];
  for (int j=0; j<16; j++)
  {
    /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
    for (int k=0; k<8; k++)
    {
      O[k] = g_aiT16[ 1][k]*src[ 16] + g_aiT16[ 3][k]*src[ 3*16] + g_aiT16[ 5][k]*src[ 5*16] + g_aiT16[ 7][k]*src[ 7*16] +
        g_aiT16[ 9][k]*src[ 9*16] + g_aiT16[11][k]*src[11*16] + g_aiT16[13][k]*src[13*16] + g_aiT16[15][k]*src[15*16];
    }
    // O_vec: the sum of every corresponding elements in odd row: O_vec[i] = g[1][i]+g[3][i]+g[5][i]...., 0<=i<16
    
    for (int k=0; k<4; k++)
    {
      EO[k] = g_aiT16[ 2][k]*src[ 2*16] + g_aiT16[ 6][k]*src[ 6*16] + g_aiT16[10][k]*src[10*16] + g_aiT16[14][k]*src[14*16];
    }
    // EO_vec: the sum of every corresponding elements in row 2, 6, 10, 14: EO_vec[i] = g[2][i]+g[6][i]+g[10][i]....,0<=i<4
    // Note: only first 4 elements will be used

    for (int k = 0; k < 2; ++k)
    {
      EEE[k] = g_aiT16[0][k]*src[ 0*16 ] + g_aiT16[8][k]*src[ 8*16 ];
      EEO[k] = g_aiT16[4][k]*src[ 4*16 ] + g_aiT16[12][k]*src[ 12*16 ];
    }
    // EEE_vec: the sum of every corresponding elements in row 0 and 8: EEE_vec[i] = g[0][i]+g[8][i], 0<=i<2
    // EEO_vec: the sum of every corresponding elements in row 4 and 12: EEO_vec[i] = g[4][i]+g[12][i], 0<=i<2

    /* Combining even and odd terms at each hierarchy levels to calculate the final spatial domain vector */
    for (int k=0; k<2; k++)
    {
      EE[k] = EEE[k] + EEO[k];
      //EE[k] = (g[0] + g[4] + g[8] + g[12])[k]; Note: only the first two elements will be used, EE[i] = g[0][i]+g[4][i], 0<=k<2
      EE[k+2] = EEE[1-k] - EEO[1-k];
      //EE[k+2] = (g[0] + g[8])- (g[4] + g[12])[1-k]; only the first two elements'sum, 0<=k<2, 
      //Note: EE[2]->EEE[1]-EEO[1], EE[3]->EEE[0]-EEO[0]
    }

    for (int k=0; k<4; k++)
    {
      E[k] = EE[k] + EO[k];
      // E[0] = EE[0] + EO[0]: the sum of the 1st elements in rows (0+8+4+12), plus the sum of the 1st elements in rows 2+6+10+14
      // E[1] = EE[1] + EO[1]: the sum of the 2nd elements in rows (0+8+4+12), plus the sum of the 2nd elements in rows 2+6+10+14
      // E[2] = EE[2] + EO[2]: the sum of the 2nd elements in rows (0+8-4-12) plus the sum of the 3rd elements in rows 2+6+10+14 
      // E[3] = EE[3] + EO[3]: the sum of the 1st elements in rows (0+8-4-12) plus the sum of 4th elements in rows 2+6+10+14 
      E[k+4] = EE[3-k] - EO[3-k];
      // E[0+4] = EE[3-0] - E0[3-0]: the sum of the 1st elements in row (0+8-4-12), minus the sum of the 4th elements in rows 2+6+10+14
      // E[1+4] = EE[3-1] - EO[3-1]: the sum of the 2nd elements in rows (0+8-4-12), minus the sum of the 3th elements in rows 2+6+10+14 
      // E[2+4] = EE[3-2] - EO[3-2]: the sum of the 2nd elements in rows (0+8+4+12), minus the sum of the 2nd elements in rows 2+6+10+14
      // E[3+4] = EE[3-3] - EO[3-3]: the sum of the 1st elements in rows (0+8+4+12), minus the sum of the 1st elements in rows 2+6+10+14
    }

    // E[0]: the 1st elements of vector(0+8+4+12+2+6+10+14) 
    // E[1]: the 2nd elements of vector(0+8+4+12+2+6+10+14) 

    // for E[2-5]: load first four elements in rows(2+6+10+14) in reverse order
    // E[2]: the 2nd element of vector (0+8-4-12)+R(2+6+10+14)
    // E[3]: the 1st element of vector (0+8-4-12)+R(2+6+10+14))
    // E[4]: the 1st element of vector (0+8-4-12)-R(2+6+10+14)
    // E[5]: the 2nd element of vector (0+8-4-12)-R(2+6+10+14)
    // for E[2-5]: load first four elements in rows(2+6+10+14) in reverse order

    // E[6]: the 2nd element of vector (0+8+4+12)-(2+6+10+14)
    // E[7]: the 1st element of vector (0+8-4-12)-(2+6+10+14)

    for (int k=0; k<8; k++)
    {
      dst[k]   = MAX( -32768, MIN( 32767, (E[k]   + O[k]   + add)>>shift ));
      dst[k+8] = MAX( -32768, MIN( 32767, (E[7-k] - O[7-k] + add)>>shift ));
      //dst[7] = (E[7]   + O[7]   + add)>>shift;
      // the 1st element of vector (0+8-4-12)-(2+6+10+14) plus the 8th element of O
      // the 1st element of vector (0+8-4-12)-(2+6+10+14) plus the 1st element of R(O)
      //dst[0+8] = (E[7-0] - O[7-0] + add)>>shift;
      // the 1st element of vector (0+8-4-12)-(2+6+10+14) minus the 1st element of R(O)

      //dst[6] = (E[6]   + O[6]   + add)>>shift;
      // the 2nd element of vector (0+8+4+12)-(2+6+10+14) plus the 7th element of O
      // the 2nd element of vector (0+8+4+12)-(2+6+10+14) plus the 2nd element of R(0)
      //dst[1+8] = (E[7-1] - O[7-1] + add)>>shift;
      // the 2nd element of vector (0+8+4+12)-(2+6+10+14) minus the 2nd element of R(0)

      //dst[5] = (E[5]   + O[5]   + add)>>shift;
      // the 2nd element of vector (0+8-4-12)-R(2+6+10+14) plus the 6th element of O
      //dst[2+8] = (E[7-2] - O[7-2] + add)>>shift;
      // the 2nd element of vector (0+8-4-12)-R(2+6+10+14) minus the 6th element of O

      //dst[4] = (E[4]   + O[4]   + add)>>shift;
      // the 1st element of vector (0+8-4-12)-R(2+6+10+14) plus the 5th element of O
      //dst[3+8] = (E[7-3] - O[7-3] + add)>>shift;
      // the 1st element of vector (0+8-4-12)-R(2+6+10+14) minus the 5th element of O

      //dst[3] = (E[3]   + O[3]   + add)>>shift;
      // the 1st element of vector (0+8-4-12)+R(2+6+10+14)) plus the 4th element of O
      //dst[4+8] = (E[7-4] - O[7-4] + add)>>shift;
      // the 1st element of vector (0+8-4-12)+R(2+6+10+14)) minus the 4th element of O

      //dst[2] = (E[2]   + O[2]   + add)>>shift;
      // the 2nd element of vector (0+8-4-12)+R(2+6+10+14) plus the 3rd element of O
      //dst[5+8] = (E[7-5] - O[7-5] + add)>>shift;
      // the 2nd element of vector (0+8-4-12)+R(2+6+10+14) minus the 3rd element of O

      //dst[1] = (E[1]   + O[1]   + add)>>shift;
      // the 2nd element of vector (1+3+5+7+9+11+13+15) + (0+8+4+12+2+6+10+14)
      //dst[6+8] = (E[7-6] - O[7-6] + add)>>shift;
      // the 2nd element of vector (1+3+5+7+9+11+13+15) - (0+8+4+12+2+6+10+14)

      //dst[0] := (E[0]   + O[0]   + add)>>shift; 
      //the 1st element of vector (1+3+5+7+9+11+13+15) + (0+8+4+12+2+6+10+14)
      //dst[7+8] = (E[7-7] - O[7-7] + add)>>shift;
      // the 1st element of vector (1+3+5+7+9+11+13+15) - (0+8+4+12+2+6+10+14)
    }
    src ++;
    dst += 16;
  }
}

//vectorization code for the inverse transform
static void partialButterflyInverse16_Intri(short *src, short *dst, int shift)
{
  // short = 2 bytes = 16 bits, 128 bits = 8 short integers
  __m128i *g_aiT16_vec = (__m128i *) g_aiT16;
  __m128i min = _mm_set1_epi32(32767);
  __m128i max = _mm_set1_epi32(-32768);
  __m128i add = _mm_set1_epi32(1<<(shift-1));
  const __m128i vm = _mm_setr_epi8(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1);

  __m128i xmm0;
  __m128i xmm1;
  //(1+3+5+7+9+11+13+15)
  __m128i xmm2 = _mm_setzero_si128();
  //(0+8+4+12+2+6+10+14)
  __m128i xmm3 = _mm_setzero_si128();// first 8 integer in a row
  //(2+6+10+14)
  __m128i xmm4 = _mm_setzero_si128();// first 8 integer in a row
  //(0+8)
  __m128i xmm5 = _mm_setzero_si128();
  //(0+8-4-12)
  __m128i xmm6 = _mm_setzero_si128();// first 8 integer in a row
  //(0+8+4+12)
  __m128i xmm7 = _mm_setzero_si128();// first 8 integer in a row
        
  // Only to compute: (1+3+5+7+9+11+13+15),
  //                  (0+8+4+12)+(2+6+10+14),
  //                  (0+8+4+12)-(2+6+10+14),
  //                  (0+8-4-12)+(2+6+10+14),
  //                  (0+8-4-12)-(2+6+10+14)
  for (int j = 0; j < 16; ++j)
  {
    for (int i = 0; i < 16; ++i)
    {
      xmm0 = _mm_load_si128(&g_aiT16_vec[2*i]);// first 8 integers in a row
      xmm1 = _mm_set1_epi16(src[16*i+j]);
      //g_aiT16[i]*src[i*16];
      xmm0 = _mm_mullo_epi16(xmm0,xmm1);
      if (i%2) //odd rows
      {
        xmm2 = _mm_adds_epi16(xmm2,xmm0);
      }
      else // even rows
      {
        xmm3 = _mm_adds_epi16(xmm3,xmm0);
        if ((i==2)||(i==6)||(i==10)||(i==14)) 
        {//(2+6+10+14)
            xmm4 = _mm_adds_epi16(xmm4,xmm0);
        }
        else
        {
          if ( (i==0) || (i==8) )
          {
            xmm5 = _mm_adds_epi16(xmm5,xmm0);
          }   
          if ( (i==4) || (i==12) ) 
          { 
            //(0+8-4-12)
            xmm6 = _mm_subs_epi16(xmm5,xmm0);
            //(0+8+4+12)
            xmm7 = _mm_adds_epi16(xmm5,xmm0);
          }
        }
      }
    }
    //xmm2: first 8 elements of odd rows'sum
    //xmm3: first 8 elements of even rows'sum
    //xmm4: first 8 elements of (2+6+10+14)
    //xmm6: first 8 elements of (0+8-4-12)
    //xmm7: first 8 elements of (0+8+4+12)
    //(0+8+4+12)-(2+6+10+14)
    xmm0 = _mm_subs_epi16(xmm7,xmm4);
    //(0+8-4-12)+(2+6+10+14)
    xmm1 = _mm_adds_epi16(xmm6,xmm4);
    //(0+8-4-12)-(2+6+10+14)
    xmm5 = _mm_subs_epi16(xmm6,xmm4);

    //re_xmm4: first 8 elements of R(2+6+10+14)
    __m128i re_xmm4 = _mm_shuffle_epi8(xmm4, vm);
    //(0+8-4-12)-R(2+6+10+14)
    __m128i xmm8 = _mm_subs_epi16(xmm6,re_xmm4);
    //(0+8-4-12)+R(2+6+10+14)
    __m128i xmm9 = _mm_adds_epi16(xmm6,re_xmm4);

    //re_0 = R(xmm2)
    __m128i re_O = _mm_shuffle_epi8(xmm2, vm);
    re_O_32 = _mm_cvtepi16_epi32(re_O);
    //dst[7] = (E[7]   + O[7]   + add)>>shift;
    // the 1st element of vector (0+8-4-12)-(2+6+10+14) plus the 1st element of R(xmm2)
    __m128i xmm5_ex = _mm_cvtepi16_epi32(xmm5);
    __m128i immediate = _mm_add_epi32(xmm5_ex,re_O_32);
    immediate = _mm_add_epi32(immediate,add);
    immediate =_mm_sra_epi32(immediate,shift);
    immediate = _mm_min_epi32(immediate,min);
    immediate = _mm_max_epi32(immediate,max);
    int result = _mm_extract_epi32(immediate, 0x0);
    dst[7] = result;

    //dst[0+8] = (E[7-0] - O[7-0] + add)>>shift;
    // the 1st element of vector (0+8-4-12)-(2+6+10+14) minus the 1st element of R(xmm2)
    immediate = _mm_sub_epi32(xmm5_ex,re_O_32);
    immediate = _mm_add_epi32(immediate,add);
    immediate = _mm_sra_epi32(immediate,shift);
    immediate = _mm_min_epi32(immediate,min);
    immediate = _mm_max_epi32(immediate,max);
    result = _mm_extract_epi32(immediate, 0x0);
    dst[8] = result;

    //dst[6] = (E[6]   + O[6]   + add)>>shift;
    // the 2nd element of vector (0+8+4+12)-(2+6+10+14) plus the 2nd element of R(xmm2)
    immediate = _mm_cvtepi16_epi32(xmm0);
    immediate = _mm_add_epi32(immediate,re_O_32);
    immediate = _mm_add_epi32(immediate,add);
    immediate =_mm_sra_epi32(immediate,shift);
    immediate = _mm_min_epi32(immediate,min);
    immediate = _mm_max_epi32(immediate,max);
    result = _mm_extract_epi32(immediate, 0x1);
    dst[6] = result;

    //dst[1+8] = (E[7-1] - O[7-1] + add)>>shift;
    // the 2nd element of vector (0+8+4+12)-(2+6+10+14) minus the 2nd element of R(xmm2)
    immediate = _mm_cvtepi16_epi32(xmm0);
    immediate = _mm_sub_epi32(immediate,re_O_32);
    immediate = _mm_add_epi32(immediate,add);
    immediate =_mm_sra_epi32(immediate,shift);
    immediate = _mm_min_epi32(immediate,min);
    immediate = _mm_max_epi32(immediate,max);
    result = _mm_extract_epi32(immediate, 0x1);
    dst[9] = result;

    int xmm8_1 = _mm_extract_epi16(xmm8, 0x1);
    int xmm2_5 = _mm_extract_epi16(xmm2, 0x5);
    //dst[5] = (E[5]   + O[5]   + add)>>shift;
    // the 2nd element of vector (0+8-4-12)-R(2+6+10+14) plus the 6th element of O
    result = (xmm8_1 + xmm2_5 + add)>>shift;
    if ( result > 32767 )
    {
      result = 32767;
    }
    if ( result < -32768 )
    {
      result = -32768;
    }
    dst[5] = result;
    //dst[2+8] = (E[7-2] - O[7-2] + add)>>shift;
    // the 2nd element of vector (0+8-4-12)-R(2+6+10+14) minus the 6th element of O
    result = (xmm8_1 - xmm2_5 + add)>>shift;
    if ( result > 32767 )
    {
      result = 32767;
    }
    if ( result < -32768 )
    {
      result = -32768;
    }
    dst[10] = result;

    int xmm8_0 = _mm_extract_epi16(xmm8, 0x0);
    int xmm2_4 = _mm_extract_epi16(xmm2, 0x4);
    //dst[4] = (E[4]   + O[4]   + add)>>shift;
    // the 1st element of vector (0+8-4-12)-R(2+6+10+14) plus the 5th element of O
    result = (xmm8_0 + xmm2_4 + add)>>shift;
    if ( result > 32767 )
    {
      result = 32767;
    }
    if ( result < -32768 )
    {
      result = -32768;
    }
    dst[4] = result;
    //dst[3+8] = (E[7-3] - O[7-3] + add)>>shift;
    // the 1st element of vector (0+8-4-12)-R(2+6+10+14) minus the 5th element of O
    result = (xmm8_0 - xmm2_4 + add)>>shift;
    if ( result > 32767 )
    {
      result = 32767;
    }
    if ( result < -32768 )
    {
      result = -32768;
    }
    dst[11] = result;

    int xmm9_0 = _mm_extract_epi16(xmm9, 0x0);
    int xmm2_3 = _mm_extract_epi16(xmm2, 0x3);
    //dst[3] = (E[3]   + O[3]   + add)>>shift;
    // the 1st element of vector (0+8-4-12)+R(2+6+10+14)) plus the 4th element of O
    result = (xmm9_0 + xmm2_3 + add)>>shift;
    if ( result > 32767 )
    {
      result = 32767;
    }
    if ( result < -32768 )
    {
      result = -32768;
    }
    dst[3] = result;
    //dst[4+8] = (E[7-4] - O[7-4] + add)>>shift;
    // the 1st element of vector (0+8-4-12)+R(2+6+10+14)) minus the 4th element of O
    result = (xmm9_0 - xmm2_3 + add)>>shift;
    if ( result > 32767 )
    {
      result = 32767;
    }
    if ( result < -32768 )
    {
      result = -32768;
    }
    dst[12] = result;

    int xmm9_1 = _mm_extract_epi16(xmm9, 0x1);
    int xmm2_2 = _mm_extract_epi16(xmm2, 0x2);
    //dst[2] = (E[2]   + O[2]   + add)>>shift;
    // the 2nd element of vector (0+8-4-12)+R(2+6+10+14) plus the 3rd element of O
    result = (xmm9_1 + xmm2_2 + add)>>shift;
    if ( result > 32767 )
    {
      result = 32767;
    }
    if ( result < -32768 )
    {
      result = -32768;
    }
    dst[2] = result;
    //dst[5+8] = (E[7-5] - O[7-5] + add)>>shift;
    // the 2nd element of vector (0+8-4-12)+R(2+6+10+14) minus the 3rd element of O
    result = (xmm9_1 - xmm2_2 + add)>>shift;
    if ( result > 32767 )
    {
      result = 32767;
    }
    if ( result < -32768 )
    {
      result = -32768;
    }
    dst[13] = result;

    //dst[0] := (E[0]   + O[0]   + add)>>shift; 
    //the 1st element of vector (1+3+5+7+9+11+13+15) + (0+8+4+12+2+6+10+14)
    //dst[1] = (E[1]   + O[1]   + add)>>shift;
    // the 2nd element of vector (1+3+5+7+9+11+13+15) + (0+8+4+12+2+6+10+14)
    __m128i xmm2_ex = _mm_cvtepi16_epi32(xmm2);
    __m128i xmm3_ex = _mm_cvtepi16_epi32(xmm3);
    immediate = _mm_add_epi32(xmm2_ex,xmm3_ex);
    immediate = _mm_add_epi32(immediate,add);
    immediate =_mm_sra_epi32(immediate,shift);
    immediate = _mm_min_epi32(immediate,min);
    immediate = _mm_max_epi32(immediate,max);
    result = _mm_extract_epi32(immediate, 0x0);
    dst[0] = result;
    result = _mm_extract_epi32(immediate, 0x1);
    dst[1] = result;

    //dst[7+8] = (E[7-7] - O[7-7] + add)>>shift;
    // the 1st element of vector (1+3+5+7+9+11+13+15) - (0+8+4+12+2+6+10+14)
    //dst[6+8] = (E[7-6] - O[7-6] + add)>>shift;
    // the 2nd element of vector (1+3+5+7+9+11+13+15) - (0+8+4+12+2+6+10+14)
    immediate = _mm_sub_epi32(xmm2_ex,xmm3_ex);
    immediate = _mm_add_epi32(immediate,add);
    immediate =_mm_sra_epi32(immediate,shift);
    immediate = _mm_min_epi32(immediate,min);
    immediate = _mm_max_epi32(immediate,max);
    result = _mm_extract_epi32(immediate, 0x0);
    dst[15] = result;
    result = _mm_extract_epi32(immediate, 0x1);
    dst[14] = result;
    dst += 16;
  }
}

static void idct16_scalar(short* pCoeff, short* pDst)
{
  short tmp[ 16*16] __attribute__((aligned(16)));
  partialButterflyInverse16(pCoeff, tmp, 7);
  partialButterflyInverse16(tmp, pDst, 12);
}

/// CURRENTLY SAME CODE AS SCALAR !!
/// REPLACE HERE WITH SSE intrinsics
static void idct16_simd(short* pCoeff, short* pDst)
{
  short tmp[ 16*16] __attribute__((aligned(16)));
  partialButterflyInverse16_Intri(pCoeff, tmp, 7);
  partialButterflyInverse16_Intri(tmp, pDst, 12);
}

int main(int argc, char **argv)
{
    //allocate memory 16-byte aligned
    short *scalar_input = (short*) memalign(16, IDCT_SIZE*IDCT_SIZE*sizeof(short));
    short *scalar_output = (short *) memalign(16, IDCT_SIZE*IDCT_SIZE*sizeof(short));

    short *simd_input = (short*) memalign(16, IDCT_SIZE*IDCT_SIZE*sizeof(short));
    short *simd_output = (short *) memalign(16, IDCT_SIZE*IDCT_SIZE*sizeof(short));

    //initialize input
    printf("input array:\n");
    for(int j=0;j<IDCT_SIZE;j++){
        for(int i=0;i<IDCT_SIZE;i++){
            short value = rand()%2 ? (rand()%32768) : -(rand()%32768) ;
            scalar_input[j*IDCT_SIZE+i] = value;
            simd_input  [j*IDCT_SIZE+i] = value;
	    printf("%d\t", value);
        }
        printf("\n");
    }

    idct16_scalar(scalar_input, scalar_output);
    idct16_simd  (simd_input  , simd_output);

    //check for correctness
    compare_results (scalar_output, simd_output, "scalar and simd");

    printf("output array:\n");
    for(int j=0;j<IDCT_SIZE;j++){
        for(int i=0;i<IDCT_SIZE;i++){
	    printf("%d\t", scalar_output[j*IDCT_SIZE+i]);
        }
        printf("\n");
    }

    //Measure the performance of each kernel
    benchmark (idct16_scalar, scalar_input, scalar_output, "scalar");
    benchmark (idct16_simd, simd_input, simd_output, "simd");

    //cleanup
    free(scalar_input);    free(scalar_output);
    free(simd_input); free(simd_output);
}

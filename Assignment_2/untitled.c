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

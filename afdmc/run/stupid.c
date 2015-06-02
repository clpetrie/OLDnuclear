#include <immintrin.h>
#include <stdio.h>
#include <math.h> 
#include <sys/time.h>


// -----------------------------------------------------------------------------
// This function returns the current time, expressed as seconds since the Epoch
// -----------------------------------------------------------------------------
double getCurrentTime(){
  struct timeval curr;
  struct timezone tz;
  gettimeofday(&curr, &tz);
  double tmp = curr.tv_sec*1000000. + curr.tv_usec;
  return tmp*1e-6;
}

// -----------------------------------------------------------------------------
// Main routine
// -----------------------------------------------------------------------------
int main() {

  int i;
  srand48(0);            // seed PRNG
  double e,s;            // timestamp variables
  float *a, *b;          // data pointers
  float *pA,*pB;         // work pointer
  __m128 rA,rB;          // variables for SSE
  __m256 rA_AVX, rB_AVX; // variables for AVX

  // define vector size 
  const int vector_size = 10000000;

  // allocate memory 
  a = (float*) _mm_malloc (vector_size*sizeof(float),32);
  b = (float*) _mm_malloc (vector_size*sizeof(float),32);

  // initialize vectors //
  for(i=0;i<vector_size;i++) {
    a[i]=fabs(drand48());
    b[i]=0.0f;
  }

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Naive implementation
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  s = getCurrentTime();
  pA = __builtin_assume_aligned(a, 16);
  pB = __builtin_assume_aligned(b, 16);
  for (i=0; i<vector_size; i++){
//    b[i] = sqrtf(sqrtf(sqrtf(a[i])));
      pB[i] = sqrtf(sqrtf(sqrtf(pA[i])));
  }
  e = getCurrentTime();
  printf("%lf ms b[42] = %lf\n",(e-s)*1000,b[42]);
//  cout << (e-s)*1000 << " ms" << ", b[42] = " << b[42] << endl;

// -----------------------------------------------------------------------------
  for(i=0;i<vector_size;i++) {
    b[i]=0.0f;
  }
// -----------------------------------------------------------------------------

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// SSE2 implementation
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pA = a; pB = b;

  s = getCurrentTime();
  for (i=0; i<vector_size; i+=4){
    rA   = _mm_load_ps(pA);
    rB   = _mm_sqrt_ps(_mm_sqrt_ps(_mm_sqrt_ps(rA)));
    _mm_store_ps(pB,rB);
    pA += 4;
    pB += 4;
  }
  e = getCurrentTime();
  printf("%lf ms b[42] = %lf\n",(e-s)*1000,b[42]);

// -----------------------------------------------------------------------------
  for(i=0;i<vector_size;i++) {
    b[i]=0.0f;
  }
// -----------------------------------------------------------------------------

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// AVX implementation
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pA = a; pB = b;

  s = getCurrentTime();
  for (i=0; i<vector_size; i+=8){
    rA_AVX   = _mm256_load_ps(pA);
    rB_AVX   = _mm256_sqrt_ps(_mm256_sqrt_ps(_mm256_sqrt_ps(rA_AVX)));
    _mm256_store_ps(pB,rB_AVX);
    pA += 8;
    pB += 8;
  }
  e = getCurrentTime();
  printf("%lf ms b[42] = %lf\n",(e-s)*1000,b[42]);

  _mm_free(a);
  _mm_free(b);

  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#define BLOCK_SIZE 32
static inline double  gtod();
static inline float* random_mat( uint32_t n );
static inline float* zero_mat( uint32_t n );

static inline void mat_mult_blocking2( float *A, float *B, float *C, const unsigned int dim);
static inline void mat_mult_blocking4( float *A, float *B, float *C, const unsigned int dim);
static inline void mat_mult_blocking8( float *A, float *B, float *C, const unsigned int dim);
static inline void mat_mult_blocking16( float *A, float *B, float *C, const unsigned int dim);
static inline void mat_mult_blocking32( float *A, float *B, float *C, const unsigned int dim);

static inline void transpose_mat( float* mat, const unsigned int dim);
static inline int compareMat(float *A, float *B, const int dim);
int main( )
{
    double t_start, t_end;
    double gflops;
    uint32_t dim  = 0;
    uint32_t step = 32;
    uint32_t max  = 2048;

    printf( "Dim;TimeBlock2;GFLOPSBlock2;TimeBlock4;GFLOPSBlock4;TimeBlock8;GFLOPSBlock8;TimeBlock16;GFLOPSBlock16;TimeBlock32;GFLOPSBlock32;\n" );

    while ( dim < max )
    {
        if ( dim < 256 )       { dim += step; }
        else if ( dim < 512 )  { dim += step*2; }
        else if ( dim < 1024 ) { dim += step*4; }
        else if ( dim < 2048 ) { dim += step*8; }
        else if ( dim < 4096 ) { dim += step*16; }
        else                   { dim += step*32; }

        float* A = random_mat( dim );
        float* B = random_mat( dim );
        float* C = zero_mat( dim );
        if ( A == NULL || B == NULL || C == NULL  )
        {
            printf( "Allocation of matrix failed.\n" );
            exit( EXIT_FAILURE );
        }
// Nothing

        t_start = gtod();
         mat_mult_blocking2(A,B,C, dim);

        t_end = gtod();
        gflops = ( ( double )2 * dim * dim * dim / 1000000000.0 ) / ( t_end - t_start );

        printf("%d;%0.4f;%0.2f;", dim, t_end - t_start, gflops );
        free( C );
        C = zero_mat( dim );

// Transpose

         t_start = gtod();
         mat_mult_blocking4(A,B,C, dim);

         t_end = gtod();
         gflops = ( ( double )2 * dim * dim * dim / 1000000000.0 ) / ( t_end - t_start );

         printf("%0.4f;%0.2f;",  t_end - t_start, gflops );
         free( C );
         C = zero_mat( dim );

// Unroll 32 Transpose*/

         t_start = gtod();
         mat_mult_blocking8(A,B,C, dim);

         t_end = gtod();
         gflops = ( ( double )2 * dim * dim * dim / 1000000000.0 ) / ( t_end - t_start );

         printf("%0.4f;%0.2f;", t_end - t_start, gflops );
         free( C );
         C = zero_mat( dim );
// Unroll 32 BLOCKING
         t_start = gtod();
         mat_mult_blocking16(A,B,C, dim);

         t_end = gtod();
         gflops = ( ( double )2 * dim * dim * dim / 1000000000.0 ) / ( t_end - t_start );

         printf("%0.4f;%0.2f;", t_end - t_start, gflops );
         free( C );
         C = zero_mat( dim );
//
// Unroll 32 BLOCKING
        t_start = gtod();
         mat_mult_blocking32(A,B,C, dim);

                  t_end = gtod();
                  gflops = ( ( double )2 * dim * dim * dim / 1000000000.0 ) / ( t_end - t_start );

                  printf("%0.4f;%0.2f;\n", t_end - t_start, gflops );
 //*/
        free( A );
        free( B );
        free( C );
    }

    printf("\n");

    return EXIT_SUCCESS;
}

static inline void mat_mult_blocking2( float *A, float *B, float *C, const unsigned int dim){
    transpose_mat(A, dim);
    uint32_t i, posXA, posYA, posYB;
    uint32_t posYB_j_dim, posXA_i, posYB_j_dim_posYA, posXA_i_dim_posYA ;
#define ADDER(a) (A[posXA_i_dim_posYA+ a] * B[posYB_j_dim_posYA + a ])
    for(posXA=0; posXA<dim; posXA +=2)
    {
        for(posYA=0; posYA<dim; posYA += 2)
        {
            for(posYB=0; posYB<dim; posYB += 2)
            {
                for(i=0; i<2; ++i)
                {
                    posXA_i = posXA + i;
                    posXA_i_dim_posYA = (posXA + i) * dim + posYA;
                    for(uint32_t j=0; j<2; ++j)
                    {
                        posYB_j_dim = (posYB + j) * dim;
                        posYB_j_dim_posYA =  (posYB + j) * dim + posYA;
                            C[ posYB_j_dim  + (posXA_i)] += ADDER(0)
                                                                + ADDER(1);//*/
                    }
                }
            }
        }
    }

    transpose_mat(A, dim);
#undef ADDER
}

static inline void mat_mult_blocking4( float *A, float *B, float *C, const unsigned int dim){
    transpose_mat(A, dim);
    uint32_t i, posXA, posYA, posYB;
    uint32_t posYB_j_dim, posXA_i, posYB_j_dim_posYA, posXA_i_dim_posYA ;
#define ADDER(a) (A[posXA_i_dim_posYA+ a] * B[posYB_j_dim_posYA + a ])
    for(posXA=0; posXA<dim; posXA +=4)
    {
        for(posYA=0; posYA<dim; posYA += 4)
        {
            for(posYB=0; posYB<dim; posYB += 4)
            {
                for(i=0; i<4; ++i)
                {
                    posXA_i = posXA + i;
                    posXA_i_dim_posYA = (posXA + i) * dim + posYA;
                    for(uint32_t j=0; j<4; ++j)
                    {
                        posYB_j_dim = (posYB + j) * dim;
                        posYB_j_dim_posYA =  (posYB + j) * dim + posYA;
                            C[ posYB_j_dim  + (posXA_i)] += ADDER(0)
                                                                + ADDER(1)
                                                                + ADDER(2)
                                                                + ADDER(3);//*/
                    }
                }
            }
        }
    }

    transpose_mat(A, dim);
#undef ADDER
}

static inline void mat_mult_blocking8( float *A, float *B, float *C, const unsigned int dim){
    transpose_mat(A, dim);
    uint32_t i, posXA, posYA, posYB;
    uint32_t posYB_j_dim, posXA_i, posYB_j_dim_posYA, posXA_i_dim_posYA ;
#define ADDER(a) (A[posXA_i_dim_posYA+ a] * B[posYB_j_dim_posYA + a ])
    for(posXA=0; posXA<dim; posXA +=8)
    {
        for(posYA=0; posYA<dim; posYA += 8)
        {
            for(posYB=0; posYB<dim; posYB += 8)
            {
#pragma Loop_Optimize Unroll
                for(i=0; i<8; ++i)
                {
                    posXA_i = posXA + i;
                    posXA_i_dim_posYA = (posXA + i) * dim + posYA;
#pragma Loop_Optimize(Unroll)
                    for(uint32_t j=0; j<8; ++j)
                    {
                        posYB_j_dim = (posYB + j) * dim;
                        posYB_j_dim_posYA =  (posYB + j) * dim + posYA;
                            C[ posYB_j_dim  + (posXA_i)] += ADDER(0)
                                                                + ADDER(1)
                                                                + ADDER(2)
                                                                + ADDER(3)
                                                                + ADDER(4)
                                                                + ADDER(5)
                                                                + ADDER(6)
                                                                + ADDER(7);//*/
                    }
                }
            }
        }
    }

    transpose_mat(A, dim);
#undef ADDER
}

static inline void mat_mult_blocking16( float *A, float *B, float *C, const unsigned int dim){
    transpose_mat(A, dim);
    uint32_t i, posXA, posYA, posYB;
    uint32_t posYB_j_dim, posXA_i, posYB_j_dim_posYA, posXA_i_dim_posYA ;
#define ADDER(a) (A[posXA_i_dim_posYA+ a] * B[posYB_j_dim_posYA + a ])
    for(posXA=0; posXA<dim; posXA +=16)
    {
        for(posYA=0; posYA<dim; posYA += 16)
        {
            for(posYB=0; posYB<dim; posYB += 16)
            {
#pragma Loop_Optimize Unroll
                for(i=0; i<16; ++i)
                {
                    posXA_i = posXA + i;
                    posXA_i_dim_posYA = (posXA + i) * dim + posYA;
#pragma Loop_Optimize(Unroll)
                    for(uint32_t j=0; j<16; ++j)
                    {
                        posYB_j_dim = (posYB + j) * dim;
                        posYB_j_dim_posYA =  (posYB + j) * dim + posYA;
                            C[ posYB_j_dim  + (posXA_i)] += ADDER(0)
                                                                + ADDER(1)
                                                                + ADDER(2)
                                                                + ADDER(3)
                                                                + ADDER(4)
                                                                + ADDER(5)
                                                                + ADDER(6)
                                                                + ADDER(7)
                                                                + ADDER(8)
                                                                + ADDER(9)
                                                                + ADDER(10)
                                                                + ADDER(11)
                                                                + ADDER(12)
                                                                + ADDER(13)
                                                                + ADDER(14)
                                                                + ADDER(15);//*/
                    }
                }
            }
        }
    }

    transpose_mat(A, dim);
#undef ADDER
}


static inline void mat_mult_blocking32( float *A, float *B, float *C, const unsigned int dim){
    transpose_mat(A, dim);
    uint32_t i, posXA, posYA, posYB;
    uint32_t posYB_j_dim, posXA_i, posYB_j_dim_posYA, posXA_i_dim_posYA ;
#define ADDER(a) (A[posXA_i_dim_posYA+ a] * B[posYB_j_dim_posYA + a ])
    for(posXA=0; posXA<dim; posXA +=32)
    {
        for(posYA=0; posYA<dim; posYA += 32)
        {
            for(posYB=0; posYB<dim; posYB += 32)
            {
#pragma Loop_Optimize Unroll
                for(i=0; i<32; ++i)
                {
                    posXA_i = posXA + i;
                    posXA_i_dim_posYA = (posXA + i) * dim + posYA;
#pragma Loop_Optimize(Unroll)
                    for(uint32_t j=0; j<32; ++j)
                    {
                        posYB_j_dim = (posYB + j) * dim;
                        posYB_j_dim_posYA =  (posYB + j) * dim + posYA;
                            C[ posYB_j_dim  + (posXA_i)] += ADDER(0)
                                                                + ADDER(1)
                                                                + ADDER(2)
                                                                + ADDER(3)
                                                                + ADDER(4)
                                                                + ADDER(5)
                                                                + ADDER(6)
                                                                + ADDER(7)
                                                                + ADDER(8)
                                                                + ADDER(9)
                                                                + ADDER(10)
                                                                + ADDER(11)
                                                                + ADDER(12)
                                                                + ADDER(13)
                                                                + ADDER(14)
                                                                + ADDER(15)
                                                                + ADDER(16)
                                                                + ADDER(17)
                                                                + ADDER(18)
                                                                + ADDER(19)
                                                                + ADDER(20)
                                                                + ADDER(21)
                                                                + ADDER(22)
                                                                + ADDER(23)
                                                                + ADDER(24)
                                                                + ADDER(25)
                                                                + ADDER(26)
                                                                + ADDER(27)
                                                                + ADDER(28)
                                                                + ADDER(29)
                                                                + ADDER(30)
                                                                + ADDER(31);//*/
                    }
                }
            }
        }
    }

    transpose_mat(A, dim);
#undef ADDER
}


static inline int compareMat(float *A, float *B, const int dim){
    const int dim2 = dim*dim;
    for(int i=0; i<dim2; ++i){
        if(abs(A[i] - B[i]) > 0.00000001){
            printf("Error in %d, %d expected %1.5f got %1.5f", i%dim, i/dim, A[i], B[i]);
            return 0;
        }
    }
    return 1;
}

/**
 * @brief Transpose a matrix
 *
*/
static inline void transpose_mat( float* mat, const unsigned int dim){
  for( uint32_t i=1; i<dim; ++i){
      for( uint32_t j=0; j<i; ++j){
          float help = mat[ i * dim +j];
          mat[ i * dim + j ] = mat[ j * dim + i];
          mat[ j * dim + i ] = help;
      }
  }
}
/** @brief Get current time stamp in seconds.
 *
 *  @return         Returns current time stamp in seconds.
 */
static inline double gtod( )
{
    struct timeval act_time;
    gettimeofday( &act_time, NULL );

    return ( double )act_time.tv_sec + ( double )act_time.tv_usec / 1000000.0;
}



/** @brief Generate randomized matrix.
 *
 *  @param dim      Dimension for the generated matrix.
 *
 *  @return         Returns a pointer to the generated matrix on success, NULL
 *                  otherwise.
 */
static inline float* random_mat( uint32_t dim )
{
    float *matrix = ( float* )malloc( sizeof( float ) * dim * dim );
    if ( matrix == NULL )
    {
        return NULL;
    }

    srand( ( unsigned ) time( NULL ) );

    for ( uint32_t i = 0; i < dim * dim; ++i)
    {
        matrix[ i ] = 1 + 1.0 / (float)(rand()% 100);
    }

  return matrix;
}


/** @brief Generate zero matrix.
 *
 *  @param dim      Dimension for the generated matrix.
 *
 *  @return         Returns a pointer to the generated matrix on success, NULL
 *                  otherwise.
 */
static inline float* zero_mat( uint32_t dim )
{
    float* matrix = ( float* )malloc( sizeof( float ) * dim * dim );
    if ( matrix == NULL )
    {
        return NULL;
    }

    for ( uint32_t i = 0; i < dim * dim; ++i)
    {
        matrix[ i ] = ( float )0.0;
    }

  return matrix;
}


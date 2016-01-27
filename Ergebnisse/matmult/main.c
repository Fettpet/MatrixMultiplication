#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#define BLOCK_SIZE 32
static inline double  gtod();
static inline double* random_mat( uint32_t n );
static inline double* zero_mat( uint32_t n );
static inline void mat_mult_non_opt( double *A, double *B, double *C, const unsigned int dim);
static inline void mat_mult_transpose( double *A, double *B, double *C, const unsigned int dim);
static inline void mat_mult_unroll_transpose( double *A, double *B, double *C, const unsigned int dim);
static inline void mat_mult_blocking( double *A, double *B, double *C, const unsigned int dim);
static inline void transpose_mat( double* mat, const unsigned int dim);
static inline int compareMat(double *A, double *B, const int dim);
int main( )
{
    double t_start, t_end;
    double gflops;
    uint32_t dim  = 0;
    uint32_t step = 32;
    uint32_t max  = 2048;

    printf( "Dim;TimeNothing;GFLOPSNothin;TimeTranspose;GFLOPSTranspose;TimeUnrollTranspose;GFLOPSTranspose;TimeBlocking;GFLOPSEBlocking;\n" );

    while ( dim < max )
    {
        if ( dim < 256 )       { dim += step; }
        else if ( dim < 512 )  { dim += step*2; }
        else if ( dim < 1024 ) { dim += step*4; }
        else if ( dim < 2048 ) { dim += step*8; }
        else if ( dim < 4096 ) { dim += step*16; }
        else                   { dim += step*32; }

        double* A = random_mat( dim );
        double* B = random_mat( dim );
        double* C = zero_mat( dim );
        if ( A == NULL || B == NULL || C == NULL  )
        {
            printf( "Allocation of matrix failed.\n" );
            exit( EXIT_FAILURE );
        }

        free( A );
        free( B );
        free( C );
    }

    printf("\n");

    return EXIT_SUCCESS;
}

/**
 * given by class
 */
static inline void mat_mult_non_opt( double *A, double *B, double *C, const unsigned int dim){
        for ( uint32_t i = 0; i < dim; i++ )
        {
            for ( uint32_t j = 0; j < dim; j++ )
            {
                for ( uint32_t k = 0; k < dim; k++ )
                {
                    // C[i][j] += A[i][k] * B[k][j]
                    C[ i * dim + j ] += A[ i * dim + k ] * B[ k * dim + j ];
                }
            }
        }
}



/**
 * Mat mult with transpose
 */
static inline void mat_mult_transpose( double *restrict A, double *B, double *C, const unsigned int dim){
    transpose_mat(B, dim);
    double *restrict b = B;
    uint32_t j, i;
    uint32_t idim, jdim;
    double temp;
    for ( j = 0; j < dim; j++ )
    {
        jdim = j * dim;
        for (  i = 0; i < dim; i++ )
        {
            idim = i * dim;
            temp = 0.0;
            for (size_t k = 0; k < dim; k++ )
            {
                    temp += A[ idim + k ] * b[ jdim + k ];
            }
            C[ idim + j ] = temp;
        }
    }
    transpose_mat(B, dim);
}



/*
static inline void mat_mult_avx2(double *A, double *B, double *C, const unsigned int dim){
    transpose_mat(B, dim);
    __m256d *vecA, *vecB;
        for ( uint32_t i = 0; i < dim; i++ )
        {
            for ( uint32_t j = 0; j < dim; j++ )
            {
                for ( uint32_t k = 0; k < dim; k+= 4 )
                {
                    vecA = &(A[i*dim + k]);
                    vecB = &(B[j*dim + k]);
                    // C[i][j] += A[i][k] * B[k][j]
                    C[ i * dim + j ] += A[ i * dim + k ] * B[ j * dim + k ];
                }
            }
        }
    transpose_mat(B, dim);
}
/**
 * Mat mult with transpose and rollout
 */
static inline void mat_mult_unroll_transpose( double *A, double *B, double *C, const unsigned int dim){
    transpose_mat(B, dim);
    uint32_t i, j, k;
    uint32_t idim, jdim, idimk, jdimk;
    for (  i = 0; i < dim; i++ )
    {
        idim = i*dim;
        for (  j = 0; j < dim; j++ )
        {
            jdim= j*dim;
                for (  k = 0; k < dim; k+=32 )
                {
                    jdimk = jdim+k;
                    idimk = idim+k;
#define ADDER(a) (A[ idimk + a] * B[ jdimk + a])
                    C[ idim + j] += ADDER(0)
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

#undef ADDER

            }
        }
    }
    transpose_mat(B, dim);
}

static inline void mat_mult_blocking( double *A, double *B, double *C, const unsigned int dim){
    transpose_mat(A, dim);
    uint32_t i, posXA, posYA, posYB;
    uint32_t posYB_j_dim, posXA_i, posYB_j_dim_posYA, posXA_i_dim_posYA ;
#define ADDER(a) (A[posXA_i_dim_posYA+ a] * B[posYB_j_dim_posYA + a ])
    for(posXA=0; posXA<dim; posXA +=BLOCK_SIZE)
    {
        for(posYA=0; posYA<dim; posYA += BLOCK_SIZE)
        {
            for(posYB=0; posYB<dim; posYB += BLOCK_SIZE)
            {
#pragma Loop_Optimize Unroll
                for(i=0; i<BLOCK_SIZE; ++i)
                {
                    posXA_i = posXA + i;
                    posXA_i_dim_posYA = (posXA + i) * dim + posYA;
#pragma Loop_Optimize(Unroll)
                    for(uint32_t j=0; j<BLOCK_SIZE; ++j)
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


static inline int compareMat(double *A, double *B, const int dim){
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
static inline void transpose_mat( double* mat, const unsigned int dim){
  for( uint32_t i=1; i<dim; ++i){
      for( uint32_t j=0; j<i; ++j){
          double help = mat[ i * dim +j];
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
static inline double* random_mat( uint32_t dim )
{
    double *matrix = ( double* )malloc( sizeof( double ) * dim * dim );
    if ( matrix == NULL )
    {
        return NULL;
    }

    srand( ( unsigned ) time( NULL ) );

    for ( uint32_t i = 0; i < dim * dim; ++i)
    {
        matrix[ i ] = 1 + 1.0 / (double)(rand()% 100);
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
static inline double* zero_mat( uint32_t dim )
{
    double* matrix = ( double* )malloc( sizeof( double ) * dim * dim );
    if ( matrix == NULL )
    {
        return NULL;
    }

    for ( uint32_t i = 0; i < dim * dim; ++i)
    {
        matrix[ i ] = ( double )0.0;
    }

  return matrix;
}


#include <stdio.h>
#include <math.h>
#include <time.h>

#include "FLAME.h"

#include "cblas.h"

/* Various constants that control what gets timed */

#define TRUE 1
#define FALSE 0

void SymMatVec_unb_var1( FLA_Obj, FLA_Obj, FLA_Obj );
void SymMatVec_unb_var4( FLA_Obj, FLA_Obj, FLA_Obj );
void SymMatVec_unb_var5( FLA_Obj, FLA_Obj, FLA_Obj );

int main(int argc, char *argv[])
{
  int n, nfirst, nlast, ninc, i, ii, jj, irep,
    nrepeats;

  double
    dtime, dtime_best, 
    diff;

  FLA_Obj
    Aobj, xobj, yobj, yold, yref;
  
  /* Initialize FLAME.  Notice, we will only use the FLAME routines to time the routines  */
  FLA_Init( );

  /* Every time trial is repeated "repeat" times and the fastest run in recorded */
  printf( "%% number of repeats:" );
  scanf( "%d", &nrepeats );
  printf( "%% %d\n", nrepeats );

  /* Timing trials for matrix sizes n=nfirst to nlast in increments 
     of ninc will be performed.  Unblocked versions are only tested to
     nlast_unb */
  printf( "%% enter nfirst, nlast, ninc:" );
  scanf( "%d%d%d", &nfirst, &nlast, &ninc );
  printf( "%% %d %d %d \n", nfirst, nlast, ninc );
  fflush( stdout );

  i = 1;
  for ( n=nfirst; n<= nlast; n+=ninc ){

    /* Allocate space for the matrix and vectors */
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Aobj );
    FLA_Obj_create( FLA_DOUBLE, n, 1, 1, n, &xobj );
    FLA_Obj_create( FLA_DOUBLE, n, 1, 1, n, &yobj );
    FLA_Obj_create( FLA_DOUBLE, n, 1, 1, n, &yold );
    FLA_Obj_create( FLA_DOUBLE, n, 1, 1, n, &yref );

    /* Generate random matrix A */
    FLA_Random_matrix( Aobj );

    /* Generate random vectors x and yold */
    FLA_Random_matrix( xobj );
    FLA_Random_matrix( yold );

    for ( irep=0; irep<nrepeats; irep++ ) {
    /* Time SymMatVec1 */
      FLA_Copy( yold, yref );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute yref = A x + y where A is symmetric stored in the
	 lower triangular part of array A, by calling SymMatVec1.  The
	 result ends up in yrefp, which we will consider to be the correct
	 result. */
      SymMatVec_unb_var1( Aobj, xobj, yref );
      
      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    printf( "data_SymMatVec1( %d, 1:2 ) = [ %d %le ];\n", i, n,
	    dtime_best );
    fflush( stdout );

    /* Time SymMatVec4 */

    for ( irep=0; irep<nrepeats; irep++ ){
      /* Copy vector yold to y */
      FLA_Copy( yold, yobj );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute yref = A x + y where A is symmetric stored in the
	 lower triangular part of array A, by calling SymMatVec2( ) */
      SymMatVec_unb_var4( Aobj, xobj, yobj );

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    diff = FLA_Max_elemwise_diff( yobj, yref );
  
    printf( "data_SymMatVec4( %d, 1:3 ) = [ %d %le %le];\n", i, n,
	    dtime_best, diff  );

    fflush( stdout );

    /* Time SymMatVec5 */

    for ( irep=0; irep<nrepeats; irep++ ){
      /* Copy vector yold to y */
      FLA_Copy( yold, yobj );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute yref = A x + y where A is symmetric stored in the lower triangular part of 
	 array A, by calling SymMatVec2( ) */
      SymMatVec_unb_var5( Aobj, xobj, yobj );

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    diff = FLA_Max_elemwise_diff( yobj, yref );
  
    printf( "data_SymMatVec5( %d, 1:3 ) = [ %d %le %le];\n", i, n,
	    dtime_best, diff  );

    fflush( stdout );

    FLA_Obj_free( &Aobj );
    FLA_Obj_free( &xobj );
    FLA_Obj_free( &yobj );
    FLA_Obj_free( &yref );
    FLA_Obj_free( &yold );

    i++;
  }
  FLA_Finalize( );

  exit( 0 );
}

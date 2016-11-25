#include <stdio.h>
#include <math.h>
#include <time.h>

#include "FLAME.h"

#include "cblas.h"

#define SIDE FLA_LEFT
#define UPLO FLA_LOWER_TRIANGULAR
#define TRANS FLA_NO_TRANSPOSE
#define DIAG  FLA_NONUNIT_DIAG

/* Various constants that control what gets timed */

#define TRUE 1
#define FALSE 0

void Trmm_xxxn_var1( FLA_Obj, FLA_Obj );

int main(int argc, char *argv[])
{
  int n, nfirst, nlast, ninc, i, ii, jj, irep,
    nrepeats;

  double
    dtime, dtime_best, 
    diff;

  FLA_Obj
    Aobj, Bobj, Bold, Bref;
  
  /* Initialize FLAME.  Notice, we will only use the FLAME routines to time the routines  */
  FLA_Init( );

  /* Every time trial is repeated "repeat" times and the fastest run in recorded */
  printf( "%% number of repeats:" );
  scanf( "%d", &nrepeats );
  printf( "%% %d\n", nrepeats );

  /* Timing trials for matrix sizes n=nfirst to nlast in increments 
     of ninc will be performed. */
  printf( "%% enter nfirst, nlast, ninc:" );
  scanf( "%d%d%d", &nfirst, &nlast, &ninc );
  printf( "%% %d %d %d \n", nfirst, nlast, ninc );
  fflush( stdout );

  i = 1;
  for ( n=nfirst; n<= nlast; n+=ninc ){

    /* Allocate space for the matrices */
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Aobj );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Bobj );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Bold );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Bref );

    /* Generate random matrices A, B, and C */
    FLA_Random_matrix( Aobj );
    FLA_Random_matrix( Bold );



    for ( irep=0; irep<nrepeats; irep++ ) {
    /* Time reference implementation (from libflame) */
      FLA_Copy( Bold, Bref );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute Bref = A B or Bref = B A, 
	 where A is triangular stored in the
	 indicated triangular part of array A, by calling FLA_Trmm.  The
	 result ends up in Bref, which we will consider to be the correct
	 result. */
      FLA_Trmm( SIDE, UPLO, TRANS, DIAG, FLA_ONE, Aobj, Bref );
      
      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    printf( "data_ref( %d, 1:2 ) = [ %d %le ];\n", i, n,
	    dtime_best );
    fflush( stdout );

    /* Time your unblocked Variant 1 */

    for ( irep=0; irep<nrepeats; irep++ ){
      /* Copy vector yold to y */
      FLA_Copy( Bold, Bobj );
    
      /* start clock */
      dtime = FLA_Clock();
 
      /* Comment out the below call and call your routine instead */
      FLA_Trmm( SIDE, UPLO, TRANS, DIAG, FLA_ONE, Aobj, Bobj );
      //      Trmm_xxxn_unb_var1( Aobj, Bobj );

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    diff = FLA_Max_elemwise_diff( Bobj, Bref );
  
    printf( "data_unb_var1( %d, 1:3 ) = [ %d %le %le];\n", i, n,
	    dtime_best, diff  );

    fflush( stdout );

    FLA_Obj_free( &Aobj );
    FLA_Obj_free( &Bobj );
    FLA_Obj_free( &Bref );
    FLA_Obj_free( &Bold );

    i++;
  }
  FLA_Finalize( );

  exit( 0 );
}

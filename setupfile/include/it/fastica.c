
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include <it/fastica.h>

#include <it/math.h>
#include <it/source.h>
#include <it/random.h>
#include <it/linalg.h>
#include <it/cplx.h>
#include <it/vec.h>
#include <it/mat.h>
#include <it/io.h>


static double
min (double x, double y)
{
  if (x >= y)
    return x;
  else
    return y;
}

static int
vec_is_negative (vec v)
{
  int cmpt, i;
  cmpt = 0;
  for (i = 0; i < vec_length (v); i++)
    {
      if (v[i] < 0)
	cmpt++;
    }

  return cmpt;
}


static vec
remmean (mat vectors, mat newVectors)
{
  int i, j;
  int m = mat_height (vectors);
  int n = mat_width (vectors);
  mat mat_int = mat_new (m, n);
  vec meanValue = vec_new (m);
  vec vec_temp_1;

  for (i = 0; i < m; i++)
    {
      vec_temp_1 = mat_get_row (vectors, i);
      meanValue[i] = vec_mean (vec_temp_1);
      vec_delete (vec_temp_1);
    }

  for (j = 0; j < n; j++)
    mat_set_col (mat_int, j, meanValue);
  mat_copy (newVectors, vectors);
  mat_sub (newVectors, mat_int);

  mat_delete (mat_int);


  return meanValue;

}

static mat
selcol (mat oldMatrix, ivec maskVector)
{
  /*Selects the columns of the matrix that marked by one in the given vector.
     The maskVector is a column vector. */

  int i;
  int j = 0;
  mat newMatrix;
  vec v;

  it_assert( ivec_length( maskVector ) == mat_width( oldMatrix ), "The mask vector and matrix are of uncompatible size."); 

  /*
  if (ivec_length (maskVector) != mat_width (oldMatrix))
    printf ("The mask vector and matrix are of uncompatible size.");
  */


  newMatrix = mat_new_zeros (mat_height (oldMatrix), ivec_sum (maskVector));
  for (i = 0; i < ivec_length (maskVector); i++)
    {
      if (maskVector[i] == 1)
	{
	  v = mat_get_col (oldMatrix, i);
	  mat_set_col (newMatrix, j, v);
	  j++;
	  vec_delete (v);
	}
    }

  return newMatrix;
}

static mat
mat_tanh (mat m)
{
  int i, j;
  mat res = mat_new (mat_height (m), mat_width (m));
  for (i = 0; i < mat_height (m); i++)
    for (j = 0; j < mat_width (m); j++)
      res[i][j] = tanh (m[i][j]);


  return res;
}

static vec
vec_tanh (vec v)
{
  int i;
  vec res = vec_new (vec_length (v));
  for (i = 0; i < vec_length (v); i++)
    res[i] = tanh (v[i]);


  return res;
}

static ivec
getSamples (int max, double percentage)
{
  int i;
  double p;
  ivec mask = ivec_new (max);
  for (i = 0; i < max; i++)
    {
      p = mt19937_rand_real1 ();
      mask[i] = (p < percentage) ? 1 : 0;
    }

  return mask;
}

static mat
mat_pow_one_by_two (mat m)
{
  vec d;
  mat D;
  mat P_inv;
  mat Dtmp, Dtmp2;
  mat P = mat_new (mat_height (m), mat_width (m));

  d = mat_eig_sym (m, P);

  vec_sqrt (d);

  P_inv = mat_new_inv (P);

  D = mat_new_diag (d);

  Dtmp = mat_new_mul (P, D);
  Dtmp2 = mat_new_mul (Dtmp, P_inv);

  vec_delete (d);
  mat_delete (P);
  mat_delete (D);
  mat_delete (Dtmp);
  mat_delete (P_inv);

  return Dtmp2;
}

static int
whitenv (mat vectors, mat E, mat D, int verbose, mat newVectors,
	 mat whiteningMatrix, mat dewhiteningMatrix)
{

  int sum_diag;
  mat D_sqrt, D_sqrt_inv, covarianceMatrix, mat_temp_1;
  vec vec_temp_1;

  vec_temp_1 = mat_get_diag (D);

  sum_diag = vec_is_negative (vec_temp_1);
  vec_delete (vec_temp_1);

  it_assert( sum_diag, "Negative eigenvalues, reduce the number of dimensions in the data." );

  /*
  if (sum_diag != 0)
    {
      printf
	("[ %d ] negative eigenvalues computed from the covariance matrix.\n",
	 sum_diag);
      printf
	("These are due to rounding errors (the correct eigenvalues are probably very small).\n");
      printf
	("To correct the situation, please reduce the number of dimensions in the data\n");
      printf
	("by using the ''lastEig'' argument in function FASTICA, or ''Reduce dim.'' button\n");
      printf ("in the graphical user interface.");

      return 1;

    }
  */


  /*Calculate the whitening and dewhitening matrices (these handle dimensionality simultaneously). */

  D_sqrt = mat_clone( D );
  mat_elem_sqrt (D_sqrt);
  D_sqrt_inv = mat_clone (D_sqrt);
  mat_inv (D_sqrt_inv);



  mat_mul_transpose_right (whiteningMatrix, D_sqrt_inv, E);
  mat_mul (dewhiteningMatrix, E, D_sqrt);


  /*Project to the eigenvectors of the covariance matrix.
     Whiten the samples and reduce dimension simultaneously. */
  if (verbose)
    fprintf (stderr, "Whitening...\n");

  mat_mul (newVectors, whiteningMatrix, vectors);


  /*Print some information to user */
  if (verbose)
    {
      mat temp = mat_new_eye (mat_height (newVectors));

      covarianceMatrix = mat_cov (newVectors);

      mat_sub (covarianceMatrix, temp);
      mat_temp_1 = mat_clone( covarianceMatrix );
      mat_elem_abs (mat_temp_1);
      fprintf (stderr, "Check: covariance differs from identity by [ %f ].\n",
	      mat_max (mat_temp_1));
      mat_delete (mat_temp_1);
      mat_delete (covarianceMatrix);
      mat_delete (temp);
    }

  mat_delete (D_sqrt);
  mat_delete (D_sqrt_inv);


  return 0;


}

static ivec
PCAmat (mat vectors, int firstEig, int lastEig, int verbose, mat D, mat E)
{

  int i, nb_eigenvalues;
  int oldDimension = mat_height (vectors);
  int maxLastEig = 0;

  double rankTolerance, lowerLimitValue, higherLimitValue, sumAll, sumUsed,
    retained;
  double sum_diag = 0;

  mat covarianceMatrix;
  vec eigenvalues, diag_D;
  ivec lowerColumns, higherColumns, selectedColumns;


  if ((lastEig < 1) || (lastEig > oldDimension))
    {
      fprintf (stderr, "Illegal value [ %d ] for parameter 'lastEig'\n", lastEig);
    }

  if ((firstEig < 1) || (firstEig > lastEig))
    {
      fprintf (stderr, "Illegal value [ %d ] for parameter 'firstEig'\n", firstEig);
    }

  /*Calculate PCA */
  /*Calculate the covariance matrix. */
  if (verbose)
    fprintf (stderr, "Calculating covariance...\n");

  covarianceMatrix = mat_cov (vectors);


  /*Calculate the eigenvalues and eigenvectors of covariance matrix. */
  eigenvalues = mat_eig_sym (covarianceMatrix, E);


  diag_D = vec_clone (eigenvalues);
  mat_diag (D, eigenvalues);
  nb_eigenvalues = vec_length (eigenvalues);

  /*The rank is determined from the eigenvalues - and not directly by
     using the function rank - because function rank uses svd, which
     in some cases gives a higher dimensionality than what can be used
     with eig later on (eig then gives negative eigenvalues). */
  rankTolerance = 1e-7; 
  for (i = 0; i < vec_length (eigenvalues); i++)
    if (eigenvalues[i] > rankTolerance)
      maxLastEig++;


  it_assert( maxLastEig, "Covariance matrix eigenvalues too small. Try resizing the data matrix." );

  /*
  if (maxLastEig == 0)
    {
      printf
	("Eigenvalues of the covariance matrix are all smaller than tolerance %f.\n",
	 rankTolerance);
      printf
	("Please make sure that your data matrix contains nonzero values.\n");
      printf
	("If the values are very small,try rescaling the data matrix.\n");
      printf ("Unable to continue, aborting.");
    }
  */

  /*Sort the eigenvalues - decending. */
  vec_sort (eigenvalues);
  vec_reverse (eigenvalues);

  /*See if the user has reduced the dimension enought */
  if (lastEig > maxLastEig)
    {
      lastEig = maxLastEig;
      if (verbose)
	fprintf(stderr, "Dimension reduced to %d due to the singularity of covariance matrix\n",lastEig - firstEig + 1);
    }
  else /*Reduce the dimensionality of the problem. */
    {
      if (verbose)
	{
	  if (oldDimension == (lastEig - firstEig + 1))
	    fprintf (stderr, "Dimension not reduced.\n");
	  else
	    fprintf (stderr, "Reducing dimension...\n");
	}
    }

  /*Drop the smaller eigenvalues */
  if (lastEig < oldDimension)
    lowerLimitValue = (eigenvalues[lastEig - 1] + eigenvalues[lastEig]) / 2;
  else
    lowerLimitValue = eigenvalues[oldDimension - 1] - 1;

  lowerColumns = ivec_new_zeros (nb_eigenvalues);
  for (i = 0; i < nb_eigenvalues; i++)
    lowerColumns[i] = (diag_D[i] > lowerLimitValue) ? 1 : 0;


  /*Drop the larger eigenvalues */
  if (firstEig > 1)
    higherLimitValue =
      (eigenvalues[firstEig - 2] + eigenvalues[firstEig - 1]) / 2;
  else
    higherLimitValue = eigenvalues[0] + 1;

  higherColumns = ivec_new_zeros (nb_eigenvalues);
  for (i = 0; i < nb_eigenvalues; i++)
    higherColumns[i] = (diag_D[i] < higherLimitValue) ? 1 : 0;


  /*Combine the results from above */
  selectedColumns = ivec_new_zeros (nb_eigenvalues);
  for (i = 0; i < nb_eigenvalues; i++)
    selectedColumns[i] = (lowerColumns[i]) & (higherColumns[i]);


  /*print some info for the user */
  if (verbose)
    fprintf (stderr, "Selected [ %d ] dimensions.\n", ivec_sum (selectedColumns));

  it_assert( ivec_sum(selectedColumns)==lastEig-firstEig+1, "Selected a wrong number of dimensions" ); 

  /*
  if (ivec_sum (selectedColumns) != (lastEig - firstEig + 1))
    printf ("Selected a wrong number of dimensions.");
  */

  if (verbose)
    {
      fprintf (stderr, "Smallest remaining (non-zero) eigenvalue [ %f ]\n", eigenvalues[lastEig - 1]);
      fprintf (stderr, "Largest remaining (non-zero) eigenvalue [ %f ]\n", eigenvalues[firstEig - 1]);
      for (i = 0; i < nb_eigenvalues; i++)
	if (selectedColumns[i] == 0)
	  sum_diag += diag_D[i];
      fprintf (stderr, "Sum of removed eigenvalues [ %f ]\n", sum_diag);

    }


  /*Some more information */
  if (verbose)
    {
      sumAll = vec_sum (eigenvalues);
      sumUsed = mat_diag_sum (D);
      retained = (sumUsed / sumAll) * 100;
      fprintf (stderr, "[ %f ] %% of (non-zero) eigenvalues retained.\n", retained);
    }

  /*free memory */

  mat_delete (covarianceMatrix);
  vec_delete (eigenvalues);
  vec_delete (diag_D);
  ivec_delete (lowerColumns);
  ivec_delete (higherColumns);

  return selectedColumns;
}

static int
fpICA (mat X, mat whiteningMatrix, mat dewhiteningMatrix, char *approach,
       int numOfIC, char *g, char *finetune, int a1, int a2, int myy,
       int stabilization, double epsilon, int maxNumIterations,
       int maxFinetune, char *initState, mat guess, double sampleSize,
       int verbose, mat A, mat W)
{

  int i, j, round;
  int vectorSize = mat_height (X);
  int numSamples = mat_width (X);
  int failureLimit;
  int stabilizationEnabled = 0;
  int gOrig = 0;
  int gFine = 0;
  int finetuningEnabled;
  int usedNlinearity;
  int stroke;
  int notFine;
  int longe;
  int numFailures, endFinetuning;
  int approachMode = 2;
  int initialStateMode = 0;
  int gabba;

  double myyK;
  double myyOrig;
  double minAbsCos, minAbsCos2;
  double Beta_double;


  mat B;
  mat BOld, BOld2;
  mat mat_temp_1, mat_temp_2, mat_temp_3;
  mat Y, Gpow3, D, Xsub, Ysub, hypTan, U, Usquared, ex, gauss, dGauss, Gskew;

  vec vec_temp_1, vec_temp_2, Beta;
  vec w, wOld, wOld2, w_minus_wOld, w_plus_wOld, w_minus_wOld2, w_plus_wOld2;
  vec EXGpow3, vec_hypTan, u, u2, vec_ex, vec_gauss, vec_dGauss, EXGskew;

  ivec mask;



  /*Checking the value for approach */
  if (!strncmp (approach, "symm", 4) || !strncmp (approach, "defl", 4))
    {
      if (!strncmp (approach, "symm", 4))
	approachMode = 1;
      if (!strncmp (approach, "defl", 4))
	approachMode = 2;

      if (verbose)
	fprintf (stderr, "Used approach [ %s ].\n", approach);
    }
  else
    {
      if(verbose)
	fprintf (stderr, "Illegal value [ %s ] for parameter: 'approach' (forcing deflation)\n", approach);
      approachMode = 2; 
    }


  /*Checking the value for numOfIC */

  it_assert( vectorSize >= numOfIC, "Must have numOfIC >= dimension!" ); 

  /*
  if (vectorSize < numOfIC)
    printf ("Must have numOfIC <= Dimension!");
  */


  /*Checking the sampleSize */
  if (sampleSize > 1)
    {
      sampleSize = 1;
      if (verbose)
	fprintf (stderr, "Warning: Setting 'sampleSize' to 1.\n");
    }
  else if (sampleSize < 1)
    {
      if ((sampleSize * numSamples) < 1000)
	{
	  sampleSize = min (1000 / numSamples, 1);

	  if (verbose)
	    fprintf(stderr, "Warning: Setting 'sampleSize' to %0.3f (%f samples).\n", sampleSize, floor (sampleSize * numSamples));
	}
    }

    if (verbose && (sampleSize < 1))
      fprintf(stderr, "Using about %2.2f%% of the samples in random order in every step.\n", sampleSize * 100);


  /*Checking the value for nonlinearity. */


    if (verbose)
      fprintf (stderr, "Used nonlinearity [ %s ].\n", g);

    if (!strncmp (g, "pow3", 4))
      gOrig = 10;
    else if (!strncmp (g, "tanh", 4))
      gOrig = 20;
    else if (!strncmp (g, "gaus", 4))
      gOrig = 30;
    else if (!strncmp (g, "skew", 4))
      gOrig = 40;
    else
      {
	gOrig = 30; 
	fprintf (stderr, "Warning: Illegal value [ %s ] for parameter: 'g' (forcing 'gaus')\n", g);
      }

    if (sampleSize != 1)
      gOrig = gOrig + 2;
    if (myy != 1)
      gOrig++;
    
    finetuningEnabled = 1;
    if (!strncmp (finetune, "pow3", 4))
      gFine = 10 + 1;
    else if (!strncmp (finetune, "tanh", 4))
      gFine = 20 + 1;
    else if (!strncmp (finetune, "gaus", 4))
      gFine = 30 + 1;
    else if (!strncmp (finetune, "skew", 4))
      gFine = 40 + 1;
    else if (!strncmp (finetune, "off", 3))
      {
	if (myy != 1)
	  gFine = gOrig;
	else
	  gFine = gOrig + 1;
	finetuningEnabled = 0;
      }
    else 
      {
	gFine = 30 + 1; 
	if ( verbose )
	  fprintf (stderr, "Warning: Illegal value [ %s ] for parameter: 'finetune' (forcing 'gaus')\n", finetune);
      }

    if (verbose && finetuningEnabled)
      fprintf (stderr, "Finetuning enabled (nonlinearity: [ %s ]).\n", finetune);

    if (stabilization)
      stabilizationEnabled = 1;
    else if (!stabilization)
      {
	if (myy != 1)
	  stabilizationEnabled = 1;
	else
	  stabilizationEnabled = 0;
      }
    else 
      {
	stabilizationEnabled = 1; 
	fprintf (stderr, "Warning: Illegal value [ %d ] for parameter: 'stabilization' (forcing stabilization)\n", stabilization);
      }

    if (verbose && stabilizationEnabled)
      fprintf (stderr, "Using stabilized algorithm.\n");



    /*Some other parameters */

    myyOrig = myy;
    /*When we start fine-tuning we'll set myy = myyK * myy */
    myyK = 0.01;
    /*How many times do we try for convergence until we give up. */
    failureLimit = 5;


    usedNlinearity = gOrig;
    stroke = 0;
    notFine = 1;
    longe = 0;
    

    /*Checking the value for initial state. */
    if (!strncmp (initState, "rand", 4))
      initialStateMode = 0;
    else if (!strncmp (initState, "guess", 4))
      {
	if ((mat_height (guess) != mat_width (whiteningMatrix)) || (mat_width (guess) != numOfIC))
	  {
	    initialStateMode = 0;
	    if (verbose)
	      fprintf(stderr, "Warning: size of initial guess is incorrect. Using random initial guess.\n");
	  }
	else
	  {
	    initialStateMode = 1;
	    if (verbose)
	      fprintf (stderr, "Using initial guess.\n");
	  }
	
      }
    else
      fprintf (stderr, "Illegal value [ %s ] for parameter: 'initState'\n", initState);

    /*********************************************************************************/
    if (verbose)
      fprintf (stderr, "Starting ICA calculation...\n");

    
    /*********************************************************************************/
    /*SYMMETRIC APPROACH */
    if (approachMode)
      {
	/* set some parameters more... */
	usedNlinearity = gOrig;
	stroke = 0;
	notFine = 1;
	longe = 0;
	B = mat_new_zeros (vectorSize, numOfIC);
	BOld = mat_new_zeros (vectorSize, numOfIC);
	BOld2 = mat_new_zeros (vectorSize, numOfIC);

	mat_zeros (A);		/*Dewhitened basis vectors. */
	mat_zeros (W);
	if (!initialStateMode)
	  {
	    /*Take random orthonormal initial vectors. */
	    for (i = 0; i < vectorSize; i++)
	      for (j = 0; j < numOfIC; j++)
		B[i][j] = it_randn ();
	    mat_gs (B);
	  }
	else if (initialStateMode)
	  {
	    /*Use the given initial vector as the initial state
	      B = whiteningMatrix * guess; */
	    mat_mul (B, whiteningMatrix, guess);
	  }
	/*
	  mat_zeros (BOld);
	  mat_zeros (BOld2);
	*/

	/*This is the actual fixed-point iteration loop. */
	for (round = 1; round <= maxNumIterations + 1; round++)
	  {
	    if (round == (maxNumIterations + 1))
	      {
		fprintf (stderr, "FastICA: No convergence after %d steps\n", maxNumIterations);
		if (B != NULL)	/* a changer */
		  {
		    /*Symmetric orthogonalization. */
		    /*B = B * real(inv(B' * B)^(1/2)); */

		    mat_temp_1 = mat_new (numOfIC, numOfIC);
		    mat_mul_transpose_left (mat_temp_1, B, B);
		    mat_inv (mat_temp_1);
		    mat_temp_2 = mat_pow_one_by_two (mat_temp_1);
		    mat_temp_3 = mat_clone (B);
		    mat_mul (B, mat_temp_3, mat_temp_2);

		    mat_delete (mat_temp_1);
		    mat_delete (mat_temp_2);
		    mat_delete (mat_temp_3);

		    /*W = B' * whiteningMatrix; */
		    mat_mul_transpose_left (W, B, whiteningMatrix);

		    /*A = dewhiteningMatrix * B; */
		    mat_mul (A, dewhiteningMatrix, B);
		  }
		else
		  {
		    /*W = [];
		      A = []; */
		    mat_void (W);
		    mat_void (A);
		  }

		return 1;
	      }

	    /*Symmetric orthogonalization.
	      B = B * real(inv(B' * B)^(1/2)); */
	    mat_temp_1 = mat_new (numOfIC, numOfIC);
	    mat_mul_transpose_left (mat_temp_1, B, B);
	    mat_inv (mat_temp_1);
	    mat_temp_2 = mat_pow_one_by_two (mat_temp_1);
	    mat_temp_3 = mat_clone (B);
	    mat_mul (B, mat_temp_3, mat_temp_2);

	    mat_delete (mat_temp_1);
	    mat_delete (mat_temp_2);
	    mat_delete (mat_temp_3);

	    /*Test for termination condition. Note that we consider opposite
	      directions here as well. */

	    /*minAbsCos = Min(abs(diag(B' * BOld))); */
	    mat_temp_1 = mat_new (numOfIC, numOfIC);
	    mat_mul_transpose_left (mat_temp_1, B, BOld);
	    vec_temp_1 = mat_get_diag (mat_temp_1);
	    mat_delete (mat_temp_1);
	    vec_abs (vec_temp_1);
	    minAbsCos = vec_min (vec_temp_1);
	    vec_delete (vec_temp_1);

	    /*minAbsCos2 = min(abs(diag(B' * BOld2))); */
	    mat_temp_1 = mat_new (numOfIC, numOfIC);
	    mat_mul_transpose_left (mat_temp_1, B, BOld2);
	    vec_temp_1 = mat_get_diag (mat_temp_1);
	    mat_delete (mat_temp_1);
	    vec_abs (vec_temp_1);
	    minAbsCos2 = vec_min (vec_temp_1);
	    vec_delete (vec_temp_1);

	    if ((1.0 - minAbsCos) < epsilon)
	      {
		if (finetuningEnabled && notFine)
		  {
		    if (verbose)
		      fprintf (stderr, "Initial convergence, fine-tuning: \n");
		    notFine = 0;
		    usedNlinearity = gFine;
		    myy = myyK * myyOrig;
		    mat_zeros (BOld);
		    mat_zeros (BOld2);
		  }
		else
		  {
		    if (verbose)
		      fprintf (stderr, "Convergence after %d steps\n", round);

		    /*Calculate the de-whitened vectors. */
		    mat_mul (A, dewhiteningMatrix, B);
		    break;
		  }
	      }

	    else if (stabilizationEnabled)
	      {
		if (!stroke && ((1 - minAbsCos2) < epsilon))
		  {
		    if (verbose)
		      fprintf (stderr, "Stroke!\n");
		    stroke = myy;
		    myy = .5 * myy;
		    if (usedNlinearity % 2 == 0)
		      usedNlinearity = usedNlinearity + 1;
		  }
		else if (stroke)
		  {
		    myy = stroke;
		    stroke = 0;
		    if ((myy==1) && ((usedNlinearity % 2) != 0))
		      usedNlinearity = usedNlinearity - 1;
		  }
		else if ((longe == 0) && (round > (maxNumIterations / 2)))
		  {
		    if (verbose)
		      fprintf (stderr, "Taking long (reducing step size)\n");
		    longe = 1;
		    myy = .5 * myy;
		    if (usedNlinearity % 2 == 0)
		      usedNlinearity = usedNlinearity + 1;
		  }
	      }

	    mat_copy (BOld2, BOld);
	    mat_copy (BOld, B);

	    /*Show the progress... */
	    if (verbose)
	      {
		if (round == 1)
		  fprintf (stderr, "Step no. %d\n", round);
		else
		  fprintf (stderr, "Step no. %d, change in value of estimate: %.3f \n", round, 1 - minAbsCos);
	      }

	    switch (usedNlinearity)
	      {
	      case 10:
		/*B = (X * (( X' * B) .^ 3)) / numSamples - 3 * B; */
		mat_temp_1 = mat_clone (B);
		mat_temp_2 = mat_new (numSamples, numOfIC);
		mat_mul_transpose_left (mat_temp_2, X, B);
		mat_temp_3 = mat_clone( mat_temp_2 );
		mat_elem_pow (mat_temp_3, 3);
		
		mat_delete (mat_temp_2);

		mat_mul (B, X, mat_temp_3);

		mat_delete (mat_temp_3);

		mat_div_by (B, numSamples);
		mat_mul_by (mat_temp_1, 3);
		mat_sub (B, mat_temp_1);

		mat_delete (mat_temp_1);

		break;


	      case 11:
		/*Y = X' * B; */
		Y = mat_new (numSamples, numOfIC);
		mat_mul_transpose_left (Y, X, B);
		
		
		/*Beta = sum(Y .* Gpow3); */
		Gpow3 = mat_clone( Y );
		mat_elem_pow (Gpow3, 3);
		mat_temp_1 = mat_clone (Y);
		mat_elem_mul (mat_temp_1, Gpow3);
		Beta = mat_cols_sum (mat_temp_1);
		mat_delete (mat_temp_1);
		
		
		/*D = diag(1 ./ (Beta - 3 * numSamples)); */
		vec_temp_1 = vec_clone (Beta);
		vec_temp_2 = vec_new_set (3 * numSamples, vec_length (Beta));
		vec_sub (vec_temp_1, vec_temp_2);
		vec_set (vec_temp_2, 1);
		vec_div (vec_temp_2, vec_temp_1);
		
		vec_delete (vec_temp_1);
		
		D = mat_new_diag (vec_temp_2);
		
		vec_delete (vec_temp_2);
		
		
		/*B = B + myy * B * (Y' * Gpow3 - diag(Beta)) * D; */
		mat_temp_2 = mat_new (numOfIC, numOfIC);
		mat_mul_transpose_left (mat_temp_2, Y, Gpow3);
		mat_temp_1 = mat_new_diag (Beta);
		mat_sub (mat_temp_2, mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		mat_temp_1 = mat_new_mul (mat_temp_2, D);
		
		mat_delete (mat_temp_2);
		
		mat_temp_2 = mat_new_mul (B, mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		mat_mul_by (mat_temp_2, myy);
		mat_add (B, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_delete (Y);
		mat_delete (Gpow3);
		vec_delete (Beta);
		mat_delete (D);
		
		break;
		
		
	      case 12:
		/*Xsub=X(:, getSamples(numSamples, sampleSize)); */
		mask = getSamples (numSamples, sampleSize);
		Xsub = selcol (X, mask);
		
		/*B = (Xsub * (( Xsub' * B) .^ 3)) / size(Xsub,2) - 3 * B; */
		mat_temp_1 = mat_new (ivec_sum (mask), numOfIC);
		mat_mul_transpose_left (mat_temp_1, Xsub, B);
		
		mat_temp_3 = mat_clone (mat_temp_1);
		mat_elem_pow (mat_temp_3, 3);
		
		mat_delete (mat_temp_1);
		
		mat_temp_2 = mat_clone (B);
		mat_mul (B, Xsub, mat_temp_3);
		
		mat_delete (mat_temp_3);
		
		mat_div_by (B, ivec_sum (mask));
		mat_mul_by (mat_temp_2, 3);
		mat_sub (B, mat_temp_2);
		
		mat_delete (mat_temp_2);
		ivec_delete (mask);
		mat_delete (Xsub);
		
		break;
		
		
		
	      case 13:
		
		/*Ysub=X(:, getSamples(numSamples, sampleSize))' * B; */
		mask = getSamples (numSamples, sampleSize);
		Ysub = mat_new (ivec_sum (mask), numOfIC);
		mat_temp_1 = selcol (X, mask);
		mat_mul_transpose_left (Ysub, mat_temp_1, B);
		
		mat_delete (mat_temp_1);
		
		/*Gpow3 = Ysub .^ 3; */
		Gpow3 = mat_clone( Ysub );
		mat_elem_pow (Gpow3, 3);
		
		/*Beta = sum(Ysub .* Gpow3); */
		mat_temp_1 = mat_clone (Ysub);
		mat_elem_mul (mat_temp_1, Gpow3);
		Beta = mat_cols_sum (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		/*D = diag(1 ./ (Beta - 3 * size(Ysub', 2))); */
		vec_temp_1 = vec_clone (Beta);
		vec_temp_2 = vec_new_set (3 * ivec_sum (mask), numOfIC);
		vec_sub (vec_temp_1, vec_temp_2);
		vec_set (vec_temp_2, 1);
		vec_div (vec_temp_2, vec_temp_1);
		
		vec_delete (vec_temp_1);
		
		D = mat_new_diag (vec_temp_2);
		
		vec_delete (vec_temp_2);
		
		
		/*B = B + myy * B * (Ysub' * Gpow3 - diag(Beta)) * D; */
		mat_temp_2 = mat_new (numOfIC, numOfIC);
		mat_mul_transpose_left (mat_temp_2, Ysub, Gpow3);
		mat_temp_1 = mat_new_diag (Beta);
		mat_sub (mat_temp_2, mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		mat_temp_1 = mat_new_mul (mat_temp_2, D);
		
		mat_delete (mat_temp_2);
		
		mat_temp_2 = mat_new_mul (B, mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		mat_mul_by (mat_temp_2, myy);
		mat_add (B, mat_temp_2);
		
		mat_delete (mat_temp_2);
		ivec_delete (mask);
		mat_delete (Ysub);
		mat_delete (Gpow3);
		vec_delete (Beta);
		mat_delete (D);
		
		break;
		
		/*tanh */
		
	      case 20:
		/*hypTan = tanh(a1 * X' * B); */
		mat_temp_1 = mat_new (numSamples, numOfIC);
		mat_mul_transpose_left (mat_temp_1, X, B);
		mat_mul_by (mat_temp_1, a1);
		hypTan = mat_tanh (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		/*B = X * hypTan / numSamples - ones(size(B,1),1) * sum(1 - hypTan .^ 2) .* B / numSamples * a1; */
		mat_temp_1 = mat_new_ones (numSamples, numOfIC);
		mat_temp_2 = mat_clone (hypTan);
		mat_elem_pow (mat_temp_2, 2);
		mat_sub (mat_temp_1, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		vec_temp_1 = mat_cols_sum (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		mat_temp_1 = mat_new (vectorSize, numOfIC);
		for (i = 0; i < vectorSize; i++)
		  mat_set_row (mat_temp_1, i, vec_temp_1);
		
		vec_delete (vec_temp_1);
		
		mat_elem_mul (mat_temp_1, B);
		mat_div_by (mat_temp_1, numSamples);
		mat_mul_by (mat_temp_1, a1);
		mat_mul (B, X, hypTan);
		mat_div_by (B, numSamples);
		mat_sub (B, mat_temp_1);
		
		mat_delete (mat_temp_1);
		mat_delete (hypTan);
		
		break;
		
	      case 21:
		
		/*Y = X' * B; */
		Y = mat_new (numSamples, numOfIC);
		mat_mul_transpose_left (Y, X, B);
		
		
		/*hypTan = tanh(a1 * Y); */
		mat_temp_1 = mat_clone (Y);
		mat_mul_by (mat_temp_1, a1);
		hypTan = mat_tanh (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		/*Beta = sum(Y .* hypTan); */
		mat_temp_1 = mat_clone (Y);
		mat_elem_mul (mat_temp_1, hypTan);
		Beta = mat_cols_sum (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		/*D = diag(1 ./ (Beta - a1 * sum(1 - hypTan .^ 2))); */
		mat_temp_1 = mat_new_ones (numSamples, numOfIC);
		mat_temp_2 = mat_clone (hypTan);
		mat_elem_pow (mat_temp_2, 2);
		mat_sub (mat_temp_1, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		vec_temp_1 = mat_cols_sum (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		vec_temp_2 = vec_clone (Beta);
		vec_mul_by (vec_temp_1, a1);
		vec_sub (vec_temp_2, vec_temp_1);
		vec_set (vec_temp_1, 1);
		vec_div (vec_temp_1, vec_temp_2);
		
		vec_delete (vec_temp_2);
		
		D = mat_new_diag (vec_temp_1);
		
		vec_delete (vec_temp_1);
		
		/*B = B + myy * B * (Y' * hypTan - diag(Beta)) * D; */
		mat_temp_1 = mat_new (numOfIC, numOfIC);
		mat_mul_transpose_left (mat_temp_1, Y, hypTan);
		mat_temp_2 = mat_new_diag (Beta);
		mat_sub (mat_temp_1, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_temp_2 = mat_new_mul (mat_temp_1, D);
		
		mat_delete (mat_temp_1);
		
		mat_temp_1 = mat_new (vectorSize, numOfIC);
		mat_mul (mat_temp_1, B, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_mul_by (mat_temp_1, myy);
		mat_add (B, mat_temp_1);
		
		mat_delete (mat_temp_1);
		mat_delete (Y);
		mat_delete (hypTan);
		vec_delete (Beta);
		mat_delete (D);
		
		break;
		
	      case 22:
		/*Xsub=X(:, getSamples(numSamples, sampleSize)); */
		mask = getSamples (numSamples, sampleSize);
		Xsub = selcol (X, mask);
		
		/*hypTan = tanh(a1 * Xsub' * B); */
		mat_temp_1 = mat_new (ivec_sum (mask), numOfIC);
		mat_mul_transpose_left (mat_temp_1, Xsub, B);
		mat_mul_by (mat_temp_1, a1);
		hypTan = mat_tanh (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		
		/*B = Xsub * hypTan / size(Xsub, 2) 
		  - ones(size(B,1),1) * sum(1 - hypTan .^ 2) .* B / size(Xsub, 2) * a1;        */
		mat_temp_1 = mat_new_ones (ivec_sum (mask), numOfIC);
		mat_temp_2 = mat_clone (hypTan);
		mat_elem_pow (mat_temp_2, 2);
		mat_sub (mat_temp_1, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		vec_temp_1 = mat_cols_sum (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		mat_temp_1 = mat_new (vectorSize, numOfIC);
		for (i = 0; i < vectorSize; i++)
		  mat_set_row (mat_temp_1, i, vec_temp_1);
		
		vec_delete (vec_temp_1);
		
		mat_elem_mul (mat_temp_1, B);
		mat_div_by (mat_temp_1, ivec_sum (mask));
		mat_mul_by (mat_temp_1, a1);
		mat_mul (B, Xsub, hypTan);
		mat_div_by (B, ivec_sum (mask));
		mat_sub (B, mat_temp_1);
		
		mat_delete (mat_temp_1);
		ivec_delete (mask);
		mat_delete (hypTan);
		
		break;

	      case 23:
		
		/*Y = X(:, getSamples(numSamples, sampleSize))' * B; */
		mask = getSamples (numSamples, sampleSize);
		Y = mat_new (ivec_sum (mask), numOfIC);
		mat_temp_1 = selcol (X, mask);
		mat_mul_transpose_left (Y, mat_temp_1, B);
		
		mat_delete (mat_temp_1);
		
		/*hypTan = tanh(a1 * Y); */
		mat_temp_1 = mat_clone (Y);
		mat_mul_by (mat_temp_1, a1);
		hypTan = mat_tanh (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		/*Beta = sum(Y .* hypTan); */
		mat_temp_1 = mat_clone (Y);
		mat_elem_mul (mat_temp_1, hypTan);
		Beta = mat_cols_sum (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		/*D = diag(1 ./ (Beta - a1 * sum(1 - hypTan .^ 2))); */
		mat_temp_1 = mat_new_ones (ivec_sum (mask), numOfIC);
		mat_temp_2 = mat_clone (hypTan);
		mat_elem_pow (mat_temp_2, 2);
		mat_sub (mat_temp_1, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		vec_temp_1 = mat_cols_sum (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		vec_temp_2 = vec_clone (Beta);
		vec_mul_by (vec_temp_1, a1);
		vec_sub (vec_temp_2, vec_temp_1);
		vec_set (vec_temp_1, 1);
		vec_div (vec_temp_1, vec_temp_2);
		
		vec_delete (vec_temp_2);
		
		D = mat_new_diag (vec_temp_1);
		
		vec_delete (vec_temp_1);
		
		
		/*B = B + myy * B * (Y' * hypTan - diag(Beta)) * D; */
		mat_temp_1 = mat_new (numOfIC, numOfIC);
		mat_mul_transpose_left (mat_temp_1, Y, hypTan);
		mat_temp_2 = mat_new_diag (Beta);
		mat_sub (mat_temp_1, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_temp_2 = mat_new_mul (mat_temp_1, D);
		
		mat_delete (mat_temp_1);
		
		mat_temp_1 = mat_new (vectorSize, numOfIC);
		mat_mul (mat_temp_1, B, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_mul_by (mat_temp_1, myy);
		mat_add (B, mat_temp_1);
		
		mat_delete (mat_temp_1);
		mat_delete (Y);
		mat_delete (hypTan);
		vec_delete (Beta);
		mat_delete (D);
		
		break;
		
		
		/*gauss */
		
	      case 30:
		/*U = X' * B; */
		U = mat_new (numSamples, numOfIC);
		mat_mul_transpose_left (U, X, B);
		
		/*Usquared=U .^ 2; */
		Usquared = mat_clone( U );
		mat_elem_pow (Usquared, 2);
		
		/*ex = exp(-a2 * Usquared / 2); */
		mat_temp_1 = mat_clone (Usquared);
		mat_mul_by (mat_temp_1, -a2);
		mat_div_by (mat_temp_1, 2);
		ex = mat_clone(mat_temp_1); 
		mat_elem_exp (ex);
		
		mat_delete (mat_temp_1);
		
		/*gauss =  U .* ex; */
		gauss = mat_clone (U);
		mat_elem_mul (gauss, ex);
		
		/*dGauss = (1 - a2 * Usquared) .*ex; */
		dGauss = mat_new_ones (numSamples, numOfIC);
		mat_temp_1 = mat_clone (Usquared);
		mat_mul_by (mat_temp_1, a2);
		mat_sub (dGauss, mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		mat_elem_mul (dGauss, ex);
		
		/*B = X * gauss / numSamples - ones(size(B,1),1) * sum(dGauss).* B / numSamples ; */
		vec_temp_1 = mat_cols_sum (dGauss);
		mat_temp_1 = mat_new (vectorSize, numOfIC);
		for (i = 0; i < vectorSize; i++)
		  mat_set_row (mat_temp_1, i, vec_temp_1);
		
		vec_delete (vec_temp_1);
		
		mat_elem_mul (mat_temp_1, B);
		mat_div_by (mat_temp_1, numSamples);
		
		mat_mul (B, X, gauss);
		mat_div_by (B, numSamples);
		mat_sub (B, mat_temp_1);
		
		mat_delete (mat_temp_1);
		mat_delete (U);
		mat_delete (Usquared);
		mat_delete (ex);
		mat_delete (gauss);
		mat_delete (dGauss);
		
		break;
		
	      case 31:
		
		/*Y = X' * B; */
		Y = mat_new (numSamples, numOfIC);
		mat_mul_transpose_left (Y, X, B);
		
		
		/*ex = exp(-a2 * (Y .^ 2) / 2); */
		mat_temp_1 = mat_clone (Y);
		mat_temp_2 = mat_clone (mat_temp_1);
		mat_elem_pow (mat_temp_2, 2);
		
		mat_delete (mat_temp_1);
		
		mat_mul_by (mat_temp_2, -a2);
		mat_div_by (mat_temp_2, 2);
		ex = mat_clone( mat_temp_2 );
		mat_elem_exp (ex);
		
		mat_delete (mat_temp_2);
		
		
		/*gauss = Y .* ex; */
		gauss = mat_clone (Y);
		mat_elem_mul (gauss, ex);
		
		
		/*Beta = sum(Y .* gauss); */
		mat_temp_1 = mat_clone (Y);
		mat_elem_mul (mat_temp_1, gauss);
		Beta = mat_cols_sum (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		
		/*D = diag(1 ./ (Beta - sum((1 - a2 * (Y .^ 2)) .* ex))); */
		mat_temp_1 = mat_new_ones (numSamples, numOfIC);
		mat_temp_2 = mat_clone (Y);
		mat_elem_pow (mat_temp_2, 2);
		mat_mul_by (mat_temp_2, a2);
		mat_sub (mat_temp_1, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_elem_mul (mat_temp_1, ex);
		vec_temp_1 = vec_clone (Beta);
		vec_temp_2 = mat_cols_sum (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		vec_sub (vec_temp_1, vec_temp_2);
		vec_ones (vec_temp_2);
		vec_div (vec_temp_2, vec_temp_1);
		
		vec_delete (vec_temp_1);
		
		D = mat_new_diag (vec_temp_2);
		
		vec_delete (vec_temp_2);
		
		
		/*B = B + myy * B * (Y' * gauss - diag(Beta)) * D; */
		mat_temp_1 = mat_new (numOfIC, numOfIC);
		mat_mul_transpose_left (mat_temp_1, Y, gauss);
		mat_temp_2 = mat_new_diag (Beta);
		mat_sub (mat_temp_1, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_temp_2 = mat_new_mul (mat_temp_1, D);
		
		mat_delete (mat_temp_1);
		
		mat_temp_1 = mat_new_mul (B, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_mul_by (mat_temp_1, myy);
		mat_add (B, mat_temp_1);
		
		mat_delete (mat_temp_1);
		mat_delete (Y);
		mat_delete (D);
		mat_delete (ex);
		mat_delete (gauss);
		vec_delete (Beta);
		
		break;
		
	      case 32:
		/*Xsub=X(:, getSamples(numSamples, sampleSize)); */
		mask = getSamples (numSamples, sampleSize);
		Xsub = selcol (X, mask);
		
		/*U = Xsub' * B; */
		U = mat_new (ivec_sum (mask), numOfIC);
		mat_mul_transpose_left (U, Xsub, B);
		
		/*Usquared=U .^ 2; */
		Usquared = mat_clone (U);
		mat_elem_pow (Usquared, 2);
		
		/*ex = exp(-a2 * Usquared / 2); */
		mat_temp_1 = mat_clone (Usquared);
		mat_mul_by (mat_temp_1, -a2);
		mat_div_by (mat_temp_1, 2);
		ex = mat_clone( mat_temp_1 );
		mat_elem_exp (ex);
		
		/*gauss =  U .* ex; */
		gauss = mat_clone (U);
		mat_elem_mul (gauss, ex);
		
		/*dGauss = (1 - a2 * Usquared) .*ex; */
		dGauss = mat_new_ones (ivec_sum (mask), numOfIC);
		mat_temp_1 = mat_clone (Usquared);
		mat_mul_by (mat_temp_1, a2);
		mat_sub (dGauss, mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		mat_elem_mul (dGauss, ex);
		
		/*B = Xsub * gauss / size(Xsub,2) - ones(size(B,1),1) * sum(dGauss).* B / size(Xsub,2); */
		vec_temp_1 = mat_cols_sum (dGauss);
		mat_temp_1 = mat_new (vectorSize, numOfIC);
		for (i = 0; i < vectorSize; i++)
		  mat_set_row (mat_temp_1, i, vec_temp_1);
		
		vec_delete (vec_temp_1);
		
		mat_elem_mul (mat_temp_1, B);
		mat_div_by (mat_temp_1, ivec_sum (mask));
		mat_mul (B, Xsub, gauss);
		mat_div_by (B, ivec_sum (mask));
		mat_sub (B, mat_temp_1);
		
		mat_delete (mat_temp_1);
		ivec_delete (mask);
		mat_delete (Xsub);
		mat_delete (U);
		mat_delete (Usquared);
		mat_delete (gauss);
		mat_delete (dGauss);
		mat_delete (ex);
		
		break;
		
	      case 33:
		
		/*Y = X(:, getSamples(numSamples, sampleSize))' * B; */
		
		mask = getSamples (numSamples, sampleSize);
		Y = mat_new (ivec_sum (mask), numOfIC);
		mat_temp_1 = selcol (X, mask);
		mat_mul_transpose_left (Y, mat_temp_1, B);
		
		mat_delete (mat_temp_1);
		
		
		/*ex = exp(-a2 * (Y .^ 2) / 2); */
		mat_temp_1 = mat_clone (Y);
		mat_temp_2 = mat_clone (mat_temp_1);
		mat_elem_pow (mat_temp_2, 2);
		
		mat_delete (mat_temp_1);
		
		mat_mul_by (mat_temp_2, -a2);
		mat_div_by (mat_temp_2, 2);
		ex = mat_clone( mat_temp_2 );
		mat_elem_exp (ex);
		
		mat_delete (mat_temp_2);
		
		
		/*gauss = Y .* ex; */
		gauss = mat_clone (Y);
		mat_elem_mul (gauss, ex);
		
		
		/*Beta = sum(Y .* gauss); */
		mat_temp_1 = mat_clone (Y);
		mat_elem_mul (mat_temp_1, gauss);
		Beta = mat_cols_sum (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		
		/*D = diag(1 ./ (Beta - sum((1 - a2 * (Y .^ 2)) .* ex))); */
		mat_temp_1 = mat_new_ones (ivec_sum (mask), numOfIC);
		mat_temp_2 = mat_clone (Y);
		mat_elem_pow (mat_temp_2, 2);
		mat_mul_by (mat_temp_2, a2);
		mat_sub (mat_temp_1, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_elem_mul (mat_temp_1, ex);
		
		vec_temp_1 = vec_clone (Beta);
		vec_temp_2 = mat_cols_sum (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		vec_sub (vec_temp_1, vec_temp_2);
		vec_ones (vec_temp_2);
		vec_div (vec_temp_2, vec_temp_1);
		
		vec_delete (vec_temp_1);
		
		D = mat_new_diag (vec_temp_2);
		
		vec_delete (vec_temp_2);
		
		/*B = B + myy * B * (Y' * gauss - diag(Beta)) * D; */
		mat_temp_1 = mat_new (numOfIC, numOfIC);
		mat_mul_transpose_left (mat_temp_1, Y, gauss);
		mat_temp_2 = mat_new_diag (Beta);
		mat_sub (mat_temp_1, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_temp_2 = mat_new_mul (mat_temp_1, D);
		
		mat_delete (mat_temp_1);
		
		mat_temp_1 = mat_new_mul (B, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_mul_by (mat_temp_1, myy);
		mat_add (B, mat_temp_1);
		
		mat_delete (mat_temp_1);
		ivec_delete (mask);
		mat_delete (Y);
		mat_delete (ex);
		mat_delete (gauss);
		vec_delete (Beta);
		mat_delete (D);
		
		break;
		
		/*skew */
		
	      case 40:
		/*B = (X * ((X' * B) .^ 2)) / numSamples; */
		mat_temp_1 = mat_new (numSamples, numOfIC);
		mat_mul_transpose_left (mat_temp_1, X, B);
		mat_temp_2 = mat_clone (mat_temp_1);
		mat_elem_pow (mat_temp_2, 2);
		
		mat_delete (mat_temp_1);
		
		mat_mul (B, X, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_div_by (B, numSamples);
		
		
		break;
		
	      case 41:
		
		/*Y = X' * B; */
		Y = mat_new (numSamples, numOfIC);
		mat_mul_transpose_left (Y, X, B);
		
		/*Gskew = Y .^ 2; */
		Gskew = mat_clone (Y);
		mat_elem_pow (Gskew, 2);
		
		/*Beta = sum(Y .* Gskew); */
		mat_temp_1 = mat_clone (Y);
		mat_elem_mul (mat_temp_1, Gskew);
		Beta = mat_cols_sum (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		/*D = diag(1 ./ (Beta)); */
		vec_temp_1 = vec_new_ones (numOfIC);
		vec_div (vec_temp_1, Beta);
		D = mat_new_diag (vec_temp_1);
		
		vec_delete (vec_temp_1);
		
		/*B = B + myy * B * (Y' * Gskew - diag(Beta)) * D; */
		mat_temp_1 = mat_new (numOfIC, numOfIC);
		mat_mul_transpose_left (mat_temp_1, Y, Gskew);
		mat_temp_2 = mat_new_diag (Beta);
		mat_sub (mat_temp_1, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_temp_2 = mat_new_mul (mat_temp_1, D);
		
		mat_delete (mat_temp_1);
		
		mat_temp_1 = mat_new (vectorSize, numOfIC);
		mat_mul (mat_temp_1, B, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_mul_by (mat_temp_1, myy);
		mat_add (B, mat_temp_1);
		
		mat_delete (mat_temp_1);
		mat_delete (Y);
		mat_delete (Gskew);
		vec_delete (Beta);
		mat_delete (D);
		
		break;
		
	      case 42:
		/*Xsub=X(:, getSamples(numSamples, sampleSize)); */
		mask = getSamples (numSamples, sampleSize);
		Xsub = selcol (X, mask);
		
		/*B = (Xsub * ((Xsub' * B) .^ 2)) / size(Xsub,2); */
		mat_temp_1 = mat_new (ivec_sum (mask), numOfIC);
		mat_mul_transpose_left (mat_temp_1, Xsub, B);
		mat_temp_2 = mat_clone (mat_temp_1);
		mat_elem_pow (mat_temp_2, 2);
		
		mat_delete (mat_temp_1);
		
		mat_mul (B, Xsub, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_div_by (B, ivec_sum (mask));
		
		mat_delete (Xsub);
		ivec_delete (mask);
		
		break;
		
	      case 43:
		
		/*Y = X(:, getSamples(numSamples, sampleSize))' * B; */
		mask = getSamples (numSamples, sampleSize);
		mat Y = mat_new (ivec_sum (mask), numOfIC);
		mat_temp_1 = selcol (X, mask);
		mat_mul_transpose_left (Y, mat_temp_1, B);
		
		mat_delete (mat_temp_1);
		
		/*Gskew = Y .^ 2; */
		Gskew = mat_clone (Y);
	        mat_elem_pow (Gskew, 2);
		
		/*Beta = sum(Y .* Gskew); */
		mat_temp_1 = mat_clone (Y);
		mat_elem_mul (mat_temp_1, Gskew);
		Beta = mat_cols_sum (mat_temp_1);
		
		mat_delete (mat_temp_1);
		
		/*D = diag(1 ./ (Beta)); */
		vec_temp_1 = vec_new_ones (numOfIC);
		vec_div (vec_temp_1, Beta);
		D = mat_new_diag (vec_temp_1);
		
		vec_delete (vec_temp_1);
		
		/*B = B + myy * B * (Y' * Gskew - diag(Beta)) * D; */
		mat_temp_1 = mat_new (numOfIC, numOfIC);
		mat_mul_transpose_left (mat_temp_1, Y, Gskew);
		mat_temp_2 = mat_new_diag (Beta);
		mat_sub (mat_temp_1, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_temp_2 = mat_new_mul (mat_temp_1, D);
		
		mat_delete (mat_temp_1);
		
		mat_temp_1 = mat_new (vectorSize, numOfIC);
		mat_mul (mat_temp_1, B, mat_temp_2);
		
		mat_delete (mat_temp_2);
		
		mat_mul_by (mat_temp_1, myy);
		mat_add (B, mat_temp_1);
		
		mat_delete (mat_temp_1);
		mat_delete (Y);
		mat_delete (Gskew);
		vec_delete (Beta);
		ivec_delete (mask);
		mat_delete (D);
		
		break;
		
	      default:
		fprintf (stderr, "Warning: Code for desired nonlinearity not found!\n");
		break;
	      }
	  }

	/*Calculate ICA filters.
	  W = B' * whiteningMatrix; */
	mat_mul_transpose_left (W, B, whiteningMatrix);
	mat_delete (B);
	mat_delete (BOld);
	mat_delete (BOld2);
      }
    
    /************************************************************************************/

    /*DEFLATION APPROACH */
    if (approachMode == 2)
      {
	
	mat_zeros (A);
	mat_zeros (W);

	B = mat_new_zeros (vectorSize, vectorSize);
	
	w = vec_new (vectorSize);
	
	wOld = vec_new (vectorSize);
	wOld2 = vec_new (vectorSize);
	w_minus_wOld = vec_new (vectorSize);
	w_plus_wOld = vec_new (vectorSize);
	w_minus_wOld2 = vec_new (vectorSize);
	w_plus_wOld2 = vec_new (vectorSize);
	/*The search for a basis vector is repeated numOfIC times. */
	round = 0;

	numFailures = 0;

	while (round < numOfIC)
	  {
	    myy = myyOrig;
	    usedNlinearity = gOrig;
	    stroke = 0;
	    notFine = 1;
	    longe = 0;
	    endFinetuning = 0;
	    
	    /*Show the progress... */
	    if (verbose)
	      fprintf (stderr, "IC %d ", round + 1);

	    /*Take a random initial vector of lenght 1 and orthogonalize it
	      with respect to the other vectors. */
	    if (!initialStateMode)
	      {
		vec_randn (w);
	      }
	    else if (initialStateMode == 1)
	      {
		vec_temp_1 = mat_get_col (guess, round);
		w = mat_vec_mul (whiteningMatrix, vec_temp_1);
		vec_delete (vec_temp_1);
	      }

	    /*w = w - B * B' * w; */
	    mat_temp_1 = mat_new_transpose (B);
	    vec_temp_1 = mat_vec_mul (mat_temp_1, w);
	    
	    mat_delete (mat_temp_1);
	    
	    vec_temp_2 = mat_vec_mul (B, vec_temp_1);
	    
	    vec_delete (vec_temp_1);
	    
	    vec_sub (w, vec_temp_2);
	    
	    vec_delete (vec_temp_2);

	    /*w = w / norm(w); */
	    vec_div_by (w, vec_norm (w, 2));
	    
	    vec_zeros (wOld);
	    vec_zeros (wOld2);
	    
	    /*This is the actual fixed-point iteration loop.
	      for(i=1;i<=maxNumIterations + 1;i++) */
	    
	    i = 1;
	    gabba = 1;
	    while (i <= maxNumIterations + gabba)
	      {
		/* Project the vector into the space orthogonal to the space
		   spanned by the earlier found basis vectors. Note that we can do
		   the projection with matrix B, since the zero entries do not
		   contribute to the projection. */
		
		/*w = w - B * B' * w; */
		mat_temp_1 = mat_new_transpose (B);
		vec_temp_1 = mat_vec_mul (mat_temp_1, w);
		
		mat_delete (mat_temp_1);
		
		vec_temp_2 = mat_vec_mul (B, vec_temp_1);
		
		vec_delete (vec_temp_1);
		
		vec_sub (w, vec_temp_2);
		
		vec_delete (vec_temp_2);
		
		/*w = w / norm(w); */
		vec_div_by (w, vec_norm (w, 2));
		
		if (notFine)
		  {
		    if (i == maxNumIterations + 1)
		      {
			if (verbose)
			  fprintf(stderr, "\nComponent number %d did not converge in %d iterations.\n", round + 1, maxNumIterations);
			
			round--;
			numFailures++;
			if (numFailures > failureLimit)
			  {
			    if (verbose)
			      fprintf(stderr, "Warning: Too many failures to converge (%d). Giving up.\n", numFailures);
			    
			    if (round == 0)
			      {
				mat_void (A);
				mat_void (W);
				/*A=[];
				  W=[]; */
			      }
			    
			    return 3;
			  }			
			break;
		      }
		  }
		else
		  {
		    if (i >= endFinetuning)
		      vec_copy (wOld, w);	/*So the algorithm will stop on the next test... */		    
		  }

		/*Show the progress... */
		if (verbose)
		  fprintf (stderr, ".");

		/*Test for termination condition. Note that the algorithm has
		  converged if the direction of w and wOld is the same, this
		  is why we test the two cases. */
		vec_copy (w_minus_wOld, w);
		vec_copy (w_plus_wOld, w);
		vec_sub (w_minus_wOld, wOld);
		vec_add (w_plus_wOld, wOld);
		
		if ((vec_norm (w_minus_wOld, 2) < epsilon) || (vec_norm (w_plus_wOld, 2) < epsilon))
		  {
		    if (finetuningEnabled && notFine)
		      {
			if (verbose)
			  fprintf (stderr, "Initial convergence, fine-tuning: ");
			notFine = 0;
			gabba = maxFinetune;
			vec_zeros (wOld);
			vec_zeros (wOld2);
			usedNlinearity = gFine;
			myy = myyK * myyOrig;

			endFinetuning = maxFinetune + i;
		      }
		    else
		      {
			numFailures = 0;
			/*Save the vector */
			/*B(:, round) = w; */
			mat_set_col (B, round, w);

			/*Calculate the de-whitened vector. */

			/*A(:,round) = dewhiteningMatrix * w; */
			vec_temp_1 = mat_vec_mul (dewhiteningMatrix, w);
			mat_set_col (A, round, vec_temp_1);

			vec_delete (vec_temp_1);
			/*Calculate ICA filter. */

			/*W(round,:) = w' * whiteningMatrix; */
			mat_temp_1 = mat_new_transpose (whiteningMatrix);
			vec_temp_1 = mat_vec_mul (mat_temp_1, w);
			
			mat_delete (mat_temp_1);
			
			mat_set_row (W, round, vec_temp_1);
			
			vec_delete (vec_temp_1);
			
			/*Show the progress... */
			if (verbose)
			  fprintf (stderr, "computed ( %d steps ) \n", i);

			break;	/*IC ready - next... */
		      }
		  }
		else if (stabilizationEnabled)
		  {
		    vec_copy (w_minus_wOld2, w);
		    vec_copy (w_plus_wOld2, w);
		    vec_sub (w_minus_wOld2, wOld2);
		    vec_add (w_plus_wOld2, wOld2);
		    
		    if (!stroke && ((vec_norm (w_minus_wOld2, 2) < epsilon) || (vec_norm (w_plus_wOld2, 2) < epsilon)))
		      {
			stroke = myy;
			if (verbose)
			  fprintf (stderr, "Stroke!");
			myy = .5 * myy;
			if (usedNlinearity % 2 == 0)
			  usedNlinearity = usedNlinearity + 1;
		      }
		    else if (stroke)
		      {
			myy = stroke;
			stroke = 0;
			if ((myy == 1) && ((usedNlinearity % 2) != 0))
			  usedNlinearity = usedNlinearity - 1;
		      }
		    else if (notFine && (longe == 0) && (i > (maxNumIterations / 2)))
		      {
			if (verbose)
			  fprintf (stderr, "Taking long (reducing step size) ");
			longe = 1;
			myy = .5 * myy;
			if ((usedNlinearity % 2) == 0)
			  usedNlinearity = usedNlinearity + 1;
		      }
		  }

		vec_copy (wOld2, wOld);
		vec_copy (wOld, w);

		switch (usedNlinearity)
		  {
		    /*pow3 */

		  case 10:
		    /*w = (X * ((X' * w) .^ 3)) / numSamples - 3 * w; */
		    mat_temp_1 = mat_new_transpose (X);
		    vec_temp_1 = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    vec_temp_2 = vec_new_pow (vec_temp_1, 3);
		    
		    vec_delete (vec_temp_1);
		    
		    vec_temp_1 = mat_vec_mul (X, vec_temp_2);
		    
		    vec_delete (vec_temp_2);
		    
		    vec_div_by (vec_temp_1, numSamples);
		    vec_mul_by (w, 3);
		    vec_sub (vec_temp_1, w);
		    vec_copy (w, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    break;
		    
		  case 11:
		    /*EXGpow3 = (X * ((X' * w) .^ 3)) / numSamples; */
		    mat_temp_1 = mat_new_transpose (X);
		    vec_temp_1 = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    vec_temp_2 = vec_new_pow (vec_temp_1, 3);
		    
		    vec_delete (vec_temp_1);
		    
		    EXGpow3 = mat_vec_mul (X, vec_temp_2);
		    
		    vec_delete (vec_temp_2);
		    
		    vec_div_by (EXGpow3, numSamples);
		    
		    /*Beta = w' * EXGpow3; */
		    Beta_double = vec_inner_product (w, EXGpow3);
		    
		    /*w = w - myy * (EXGpow3 - Beta * w) / (3 - Beta); */
		    vec_temp_1 = vec_clone (w);
		    vec_mul_by (vec_temp_1, Beta_double);
		    vec_sub (EXGpow3, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    vec_mul_by (EXGpow3, myy);
		    vec_div_by (EXGpow3, (3.0 - Beta_double));
		    vec_sub (w, EXGpow3);
		    
		    vec_delete (EXGpow3);
		    
		    break;
		    
		  case 12:
		    /*Xsub=X(:,getSamples(numSamples, sampleSize)); */
		    mask = getSamples (numSamples, sampleSize);
		    Xsub = selcol (X, mask);
		    
		    /*w = (Xsub * ((Xsub' * w) .^ 3)) / size(Xsub, 2) - 3 * w; */
		    mat_temp_1 = mat_new_transpose (Xsub);
		    vec_temp_1 = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    vec_temp_2 = vec_new_pow (vec_temp_1, 3);
		    
		    vec_delete (vec_temp_1);
		    
		    vec_temp_1 = mat_vec_mul (Xsub, vec_temp_2);
		    
		    vec_delete (vec_temp_2);
		    
		    vec_div_by (vec_temp_1, ivec_sum (mask));
		    vec_mul_by (w, 3);
		    vec_sub (vec_temp_1, w);
		    vec_copy (w, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    ivec_delete (mask);
		    mat_delete (Xsub);
		    
		    break;
		    
		  case 13:
		    /*Xsub=X(:,getSamples(numSamples, sampleSize)); */
		    mask = getSamples (numSamples, sampleSize);
		    Xsub = selcol (X, mask);
		    
		    /*EXGpow3 = (Xsub * ((Xsub' * w) .^ 3)) / size(Xsub, 2); */
		    mat_temp_1 = mat_new_transpose (Xsub);
		    vec_temp_1 = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    vec_temp_2 = vec_new_pow (vec_temp_1, 3);
		    
		    vec_delete (vec_temp_1);
		    
		    EXGpow3 = mat_vec_mul (Xsub, vec_temp_2);
		    
		    vec_delete (vec_temp_2);
		    
		    vec_div_by (EXGpow3, ivec_sum (mask));
		    
		    /*Beta = w' * EXGpow3; */
		    Beta_double = vec_inner_product (w, EXGpow3);
		    
		    /*w = w - myy * (EXGpow3 - Beta * w) / (3 - Beta); */
		    vec_temp_1 = vec_clone (w);
		    vec_mul_by (vec_temp_1, Beta_double);
		    
		    vec_sub (EXGpow3, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    vec_mul_by (EXGpow3, myy);
		    vec_div_by (EXGpow3, (3.0 - Beta_double));
		    vec_sub (w, EXGpow3);
		    
		    mat_delete (Xsub);
		    ivec_delete (mask);
		    vec_delete (EXGpow3);
		    
		    break;
		    
		    /*tanh */
		    
		  case 20:
		    /*hypTan = tanh(a1 * X' * w); */
		    mat_temp_1 = mat_new_transpose (X);
		    vec_temp_1 = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    vec_mul_by (vec_temp_1, a1);
		    vec_hypTan = vec_tanh (vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    /*w = (X * hypTan - a1 * sum(1 - hypTan .^ 2)' * w) / numSamples; */
		    vec_temp_1 = mat_vec_mul (X, vec_hypTan);
		    
		    vec_mul_by (w, numSamples - vec_sum_sqr (vec_hypTan));
		    vec_mul_by (w, a1);
		    
		    vec_sub (vec_temp_1, w);
		    vec_div_by (vec_temp_1, numSamples);
		    vec_copy (w, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    vec_delete (vec_hypTan);
		    
		    break;
		    
		  case 21:
		    /*hypTan = tanh(a1 * X' * w); */
		    mat_temp_1 = mat_new_transpose (X);
		    vec_temp_1 = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    vec_mul_by (vec_temp_1, a1);
		    vec_hypTan = vec_tanh (vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    /*Beta = w' * X * hypTan; */
		    vec_temp_1 = mat_vec_mul (X, vec_hypTan);
		    Beta_double = vec_inner_product (w, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    /*w = w - myy * ((X * hypTan - Beta * w) / (a1 * sum((1-hypTan .^2)') - Beta)); */
		    vec_temp_1 = vec_clone (w);
		    vec_mul_by (vec_temp_1, Beta_double);
		    vec_temp_2 = mat_vec_mul (X, vec_hypTan);
		    vec_sub (vec_temp_2, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    vec_mul_by (vec_temp_2, myy);
		    vec_div_by (vec_temp_2,
				(a1 * (numSamples - vec_sum_sqr (vec_hypTan)) -
				 Beta_double));
		    vec_sub (w, vec_temp_2);
		    
		    vec_delete (vec_temp_2);
		    vec_delete (vec_hypTan);

		    break;
		    
		  case 22:
		    /*Xsub=X(:,getSamples(numSamples, sampleSize)); */
		    mask = getSamples (numSamples, sampleSize);
		    Xsub = selcol (X, mask);
		    
		    /*hypTan = tanh(a1 * Xsub' * w); */
		    mat_temp_1 = mat_new_transpose (Xsub);
		    vec_temp_1 = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    vec_mul_by (vec_temp_1, a1);
		    vec_hypTan = vec_tanh (vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    /*w = (Xsub * hypTan - a1 * sum(1 - hypTan .^ 2)' * w) / size(Xsub, 2); */
		    vec_temp_1 = mat_vec_mul (Xsub, vec_hypTan);
		    
		    
		    vec_mul_by (w, ivec_sum (mask) - vec_sum_sqr (vec_hypTan));
		    vec_mul_by (w, a1);
		    
		    vec_sub (vec_temp_1, w);
		    vec_div_by (vec_temp_1, ivec_sum (mask));
		    vec_copy (w, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    ivec_delete (mask);
		    mat_delete (Xsub);
		    vec_delete (vec_hypTan);
		    
		    break;
		    
		  case 23:
		    /*Xsub=X(:,getSamples(numSamples, sampleSize)); */
		    mask = getSamples (numSamples, sampleSize);
		    Xsub = selcol (X, mask);
		    
		    /*hypTan = tanh(a1 * Xsub' * w); */
		    mat_temp_1 = mat_new_transpose (Xsub);
		    vec_temp_1 = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    vec_mul_by (vec_temp_1, a1);
		    vec_hypTan = vec_tanh (vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    /*Beta = w' * Xsub * hypTan; */
		    vec_temp_1 = mat_vec_mul (Xsub, vec_hypTan);
		    Beta_double = vec_inner_product (w, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    /*w = w - myy * ((Xsub * hypTan - Beta * w) / (a1 * sum((1-hypTan .^2)') - Beta)); */
		    vec_temp_1 = vec_clone (w);
		    vec_mul_by (vec_temp_1, Beta_double);
		    vec_temp_2 = mat_vec_mul (Xsub, vec_hypTan);
		    vec_sub (vec_temp_2, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    vec_mul_by (vec_temp_2, myy);
		    vec_div_by (vec_temp_2,
				(a1 *
				 (ivec_sum (mask) - vec_sum_sqr (vec_hypTan)) -
				 Beta_double));
		    vec_sub (w, vec_temp_2);
		    
		    vec_delete (vec_temp_2);
		    
		    mat_delete (Xsub);
		    ivec_delete (mask);
		    vec_delete (vec_hypTan);
		    
		    break;
		    
		    /*gauss */
		    
		  case 30:
		    //This has been split for performance reasons.
		    /*u = X' * w; */
		    mat_temp_1 = mat_new_transpose (X);
		    u = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    /*u2=u.^2; */
		    u2 = vec_new_pow (u, 2);
		    
		    /*ex=exp(-a2 * u2/2); */
		    vec_ex = vec_clone (u2);
		    vec_mul_by (vec_ex, -a2);
		    vec_div_by (vec_ex, 2);
		    vec_exp (vec_ex);
		    
		    
		    /*gauss =  u.*ex; */
		    vec_gauss = vec_clone (u);
		    vec_mul (vec_gauss, vec_ex);
		    
		    /*dGauss = (1 - a2 * u2) .*ex; */
		    vec_dGauss = vec_new_ones (numSamples);
		    vec_mul_by (u2, a2);
		    vec_sub (vec_dGauss, u2);
		    vec_mul (vec_dGauss, vec_ex);
		    
		    /*w = (X * gauss - sum(dGauss)' * w) / numSamples; */
		    vec_temp_1 = mat_vec_mul (X, vec_gauss);
		    vec_mul_by (w, vec_sum (vec_dGauss));
		    vec_sub (vec_temp_1, w);
		    vec_div_by (vec_temp_1, numSamples);
		    vec_copy (w, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    vec_delete (u);
		    vec_delete (u2);
		    vec_delete (vec_ex);
		    vec_delete (vec_gauss);
		    vec_delete (vec_dGauss);
		    
		    break;
		    
		  case 31:
		    /*u = X' * w; */
		    mat_temp_1 = mat_new_transpose (X);
		    u = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    /*u2=u.^2; */
		    u2 = vec_new_pow (u, 2);
		    
		    /*ex=exp(-a2 * u2/2); */
		    vec_ex = vec_clone (u2);
		    vec_mul_by (vec_ex, -a2);
		    vec_div_by (vec_ex, 2);
		    vec_exp (vec_ex);
		    
		    /*gauss =  u.*ex; */
		    vec_gauss = vec_clone (u);
		    vec_mul (vec_gauss, vec_ex);
		    
		    /*dGauss = (1 - a2 * u2) .*ex; */
		    vec_dGauss = vec_new_ones (numSamples);
		    vec_mul_by (u2, a2);
		    vec_sub (vec_dGauss, u2);
		    vec_mul (vec_dGauss, vec_ex);
		    
		    /*Beta = w' * X * gauss; */
		    vec_temp_1 = mat_vec_mul (X, vec_gauss);
		    Beta_double = vec_inner_product (w, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    /*w = w - myy * ((X * gauss - Beta * w) / (sum(dGauss)' - Beta)); */
		    vec_temp_1 = vec_clone (w);
		    vec_mul_by (vec_temp_1, Beta_double);
		    vec_temp_2 = mat_vec_mul (X, vec_gauss);
		    vec_sub (vec_temp_2, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    vec_mul_by (vec_temp_2, myy);
		    vec_div_by (vec_temp_2,
				(vec_sum (vec_dGauss) - Beta_double));
		    vec_sub (w, vec_temp_2);
		    
		    vec_delete (vec_temp_2);
		    vec_delete (u);
		    vec_delete (u2);
		    vec_delete (vec_ex);
		    vec_delete (vec_gauss);
		    vec_delete (vec_dGauss);
		    
		    break;
		    
		  case 32:
		    /*Xsub=X(:,getSamples(numSamples, sampleSize)); */
		    mask = getSamples (numSamples, sampleSize);
		    Xsub = selcol (X, mask);
		    
		    /*u = Xsub' * w; */
		    mat_temp_1 = mat_new_transpose (Xsub);
		    u = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    /*u2=u.^2; */
		    u2 = vec_new_pow (u, 2);
		    
		    /*ex=exp(-a2 * u2/2); */
		    vec_ex = vec_clone (u2);
		    vec_mul_by (vec_ex, -a2);
		    vec_div_by (vec_ex, 2);
		    vec_exp (vec_ex);
		    
		    /*gauss =  u.*ex; */
		    vec_gauss = vec_clone (u);
		    vec_mul (vec_gauss, vec_ex);
		    
		    /*dGauss = (1 - a2 * u2) .*ex; */
		    vec_dGauss = vec_new_ones (ivec_sum (mask));
		    vec_mul_by (u2, a2);
		    vec_sub (vec_dGauss, u2);
		    vec_mul (vec_dGauss, vec_ex);
		    
		    /*w = (Xsub * gauss - sum(dGauss)' * w) / size(Xsub, 2); */
		    vec_temp_1 = mat_vec_mul (Xsub, vec_gauss);
		    vec_mul_by (w, vec_sum (vec_dGauss));
		    vec_sub (vec_temp_1, w);
		    vec_div_by (vec_temp_1, ivec_sum (mask));
		    vec_copy (w, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    ivec_delete (mask);
		    mat_delete (Xsub);
		    vec_delete (u);
		    vec_delete (u2);
		    vec_delete (vec_ex);
		    vec_delete (vec_gauss);
		    vec_delete (vec_dGauss);
		    
		    break;
		    
		  case 33:
		    /*Xsub=X(:,getSamples(numSamples, sampleSize)); */
		    mask = getSamples (numSamples, sampleSize);
		    Xsub = selcol (X, mask);
		    
		    /*u = Xsub' * w; */
		    mat_temp_1 = mat_new_transpose (Xsub);
		    u = mat_vec_mul (mat_temp_1, w);
		    mat_delete (mat_temp_1);
		    
		    /*u2=u.^2; */
		    u2 = vec_new_pow (u, 2);
		    
		    /*ex=exp(-a2 * u2/2); */
		    vec_ex = vec_clone (u2);
		    vec_mul_by (vec_ex, -a2);
		    vec_div_by (vec_ex, 2);
		    vec_exp (vec_ex);
		    
		    /*gauss =  u.*ex; */
		    vec_gauss = vec_clone (u);
		    vec_mul (vec_gauss, vec_ex);
		    
		    /*dGauss = (1 - a2 * u2) .*ex; */
		    vec_dGauss = vec_new_ones (ivec_sum (mask));
		    vec_mul_by (u2, a2);
		    vec_sub (vec_dGauss, u2);
		    vec_mul (vec_dGauss, vec_ex);
		    
		    /*Beta = w' * Xsub * gauss; */
		    vec_temp_1 = mat_vec_mul (Xsub, vec_gauss);
		    Beta_double = vec_inner_product (w, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    /*w = w - myy * ((Xsub * gauss - Beta * w) / (sum(dGauss)' - Beta)); */
		    vec_temp_1 = vec_clone (w);
		    vec_mul_by (vec_temp_1, Beta_double);
		    vec_temp_2 = mat_vec_mul (Xsub, vec_gauss);
		    vec_sub (vec_temp_2, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    vec_mul_by (vec_temp_2, myy);
		    vec_div_by (vec_temp_2,
				(vec_sum (vec_dGauss) - Beta_double));
		    vec_sub (w, vec_temp_2);
		    
		    vec_delete (vec_temp_2);
		    ivec_delete (mask);
		    mat_delete (Xsub);
		    vec_delete (u);
		    vec_delete (u2);
		    vec_delete (vec_ex);
		    vec_delete (vec_gauss);
		    vec_delete (vec_dGauss);
		    
		    break;
		    
		    /*skew */
		    
		  case 40:
		    /*w = (X * ((X' * w) .^ 2)) / numSamples; */
		    mat_temp_1 = mat_new_transpose (X);
		    vec_temp_1 = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    vec_temp_2 = vec_new_pow (vec_temp_1, 2);
		    
		    vec_delete (vec_temp_1);
		    
		    vec_temp_1 = mat_vec_mul (X, vec_temp_2);
		    vec_copy (w, vec_temp_2);
		    
		    vec_delete (vec_temp_2);
		    
		    vec_div_by (w, numSamples);
		    
		    break;
		    
		  case 41:
		    /*EXGskew = (X * ((X' * w) .^ 2)) / numSamples; */
		    mat_temp_1 = mat_new_transpose (X);
		    vec_temp_1 = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    vec_temp_2 = vec_new_pow (vec_temp_1, 2);
		    
		    vec_delete (vec_temp_1);
		    
		    EXGskew = mat_vec_mul (X, vec_temp_2);
		    
		    vec_delete (vec_temp_2);
		    
		    vec_div_by (EXGskew, numSamples);
		    
		    /*Beta = w' * EXGskew; */
		    Beta_double = vec_inner_product (w, EXGskew);
		    
		    /*w = w - myy * (EXGskew - Beta*w)/(-Beta); */
		    vec_temp_1 = vec_clone (w);
		    vec_mul_by (vec_temp_1, Beta_double);
		    vec_sub (EXGskew, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    vec_mul_by (EXGskew, myy);
		    vec_div_by (EXGskew, Beta_double);
		    vec_add (w, EXGskew);
		    
		    vec_delete (EXGskew);
		    
		    break;
		    
		  case 42:
		    /*Xsub=X(:,getSamples(numSamples, sampleSize)); */
		    mask = getSamples (numSamples, sampleSize);
		    Xsub = selcol (X, mask);
		    
		    /*w = (Xsub * ((Xsub' * w) .^ 2)) / size(Xsub, 2); */
		    mat_temp_1 = mat_new_transpose (Xsub);
		    vec_temp_1 = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    vec_temp_2 = vec_new_pow (vec_temp_1, 2);
		    
		    vec_delete (vec_temp_1);
		    
		    vec_temp_1 = mat_vec_mul (Xsub, vec_temp_2);
		    
		    vec_delete (vec_temp_2);
		    
		    vec_div_by (vec_temp_1, ivec_sum (mask));
		    vec_copy (w, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    ivec_delete (mask);
		    mat_delete (Xsub);
		    
		    break;
		    
		  case 43:
		    /*Xsub=X(:,getSamples(numSamples, sampleSize)); */
		    mask = getSamples (numSamples, sampleSize);
		    Xsub = selcol (X, mask);
		    
		    /*EXGskew = (Xsub * ((Xsub' * w) .^ 2)) / size(Xsub, 2); */
		    mat_temp_1 = mat_new_transpose (Xsub);
		    vec_temp_1 = mat_vec_mul (mat_temp_1, w);
		    
		    mat_delete (mat_temp_1);
		    
		    vec_temp_2 = vec_new_pow (vec_temp_1, 2);
		    
		    vec_delete (vec_temp_1);
		    
		    EXGskew = mat_vec_mul (Xsub, vec_temp_2);
		    
		    vec_delete (vec_temp_2);
		    
		    vec_div_by (EXGskew, ivec_sum (mask));
		    
		    /*Beta = w' * EXGskew; */
		    Beta_double = vec_inner_product (w, EXGskew);
		    
		    /*w = w - myy * (EXGskew - Beta*w)/(-Beta); */
		    vec_temp_1 = vec_clone (w);
		    vec_mul_by (vec_temp_1, Beta_double);
		    vec_sub (EXGskew, vec_temp_1);
		    
		    vec_delete (vec_temp_1);
		    
		    vec_mul_by (EXGskew, myy);
		    vec_div_by (EXGskew, Beta_double);
		    vec_add (w, EXGskew);

		    mat_delete (Xsub);
		    ivec_delete (mask);
		    vec_delete (EXGskew);
		    
		    break;
		    
		  case '?':
		    printf ("Code for desired nonlinearity not found!\n");
		    break;
		  }
		
		/*Normalize the new w. */
		/*w = w / norm(w); */
		vec_div_by (w, vec_norm (w, 2));
		i++;
	      }
	    
	    round++;
	  }
	
	if (verbose)
	  fprintf (stderr, "Done.\n");

	mat_delete (B);
	vec_delete (w);
	vec_delete (wOld);
	vec_delete (wOld2);
	vec_delete (w_plus_wOld);
	vec_delete (w_minus_wOld);
	vec_delete (w_plus_wOld2);
	vec_delete (w_minus_wOld2);
	
      }


    return 0;
    
}

int
it_ica (mat X, char *approach, int numOfIC, char *g, char *finetune, int a1,
	 int a2, double myy, int stabilization, double epsilon,
	 int maxNumIterations, int maxFinetune, char *initState, mat guess,
	 int sampleSize, int verbose, int jumpPCA, int jumpWhitening,
	 int only, int userNumOfIC, int firstEIG, int lastEIG, mat * icasig,
	 mat * newE, mat * newD, mat * whitesig, mat * whiteningMatrix,
	 mat * dewhiteningMatrix, mat * A, mat * W)
{

  int i;
  int Dim = mat_height (X);
  int NumOfSampl = mat_width (X);
  int newDim = 0;

  ivec selectedColumns;
  mat mat_temp_1;
  vec vec_temp_1;
  mat D = mat_new (Dim, Dim);
  mat E = mat_new (Dim, Dim);

  vec mixedmean;

  mat mixedsig = mat_new (Dim, NumOfSampl);
  mixedmean = remmean (X, mixedsig);

  if (verbose)
    {
      fprintf (stderr, "Number of signals: %d\n", mat_height (mixedsig));
      fprintf (stderr, "Number of samples: %d\n", mat_width (mixedsig));
    }

  /*Check if the data has been entered the wrong way,
    but warn only... it may be on purpose */

  if (Dim > NumOfSampl)
    {
      if (verbose)
	{
	  fprintf (stderr, "Warning: ");
	  fprintf (stderr, "The signal matrix may be oriented in the wrong way.\n");
	  fprintf (stderr, "In that case transpose the matrix.\n\n");
	}
      
      return 1;
    }

  /*PCA*/
  /*We need the results of PCA for whitening, but if we don't
    need to do whitening... then we dont need PCA... */
  if (jumpWhitening == 3)
    {
      if (verbose)
	{
	  fprintf (stderr, "Whitened signal and corresponding matrises supplied\n");
	  fprintf (stderr, "PCA calculations not needed\n");
	}
    }

  /*OK, so first we need to calculate PCA
    Check to see if we already have the PCA data */
  else
    {
      if (jumpPCA == 2)
	{
	  if (verbose)
	    {
	      fprintf (stderr, "Values for PCA calculations supplied\n");
	      fprintf (stderr, "PCA calculations not needed\n");
	    }
	}
      else
	{
	  /*display notice if the user entered one, but not both, of E and D. */
	  if ((jumpPCA > 0) && verbose)
	    {
	      fprintf (stderr, "You must suply all of these in order to jump PCA:\n");
	      fprintf (stderr, "'pcaE', 'pcaD'.\n");
	    }

	  /*Calculate PCA */
	  selectedColumns =
	    PCAmat (mixedsig, firstEIG, lastEIG, verbose, D, E);
	  newDim = ivec_sum (selectedColumns);
	  *(newE) = selcol (E, selectedColumns);
	  mat_temp_1 = selcol (D, selectedColumns);
	  mat_transpose (mat_temp_1);
	  *(newD) = selcol (mat_temp_1, selectedColumns);
	  mat_delete (mat_temp_1);
	  ivec_delete (selectedColumns);
	}
    }
  
  /*skip the rest if user only wanted PCA */
  if (only > 1)
    {      
      /*Whitening the data */
      *(whitesig) = mat_new (newDim, NumOfSampl);
      *(whiteningMatrix) = mat_new (newDim, Dim);
      *(dewhiteningMatrix) = mat_new (Dim, newDim);

      /*Check to see if the whitening is needed... */
      if (jumpWhitening == 3)
	{
	  if (verbose)
	    fprintf (stderr, "Whitening not needed.\n");
	}
      else
	{
	  /*Whitening is needed
	    display notice if the user entered some of the whitening info, but not all. */
	  if ((jumpWhitening > 0) && verbose)
	    {
	      fprintf(stderr, "You must suply all of these in order to jump whitening:\n");
	      fprintf(stderr, "'whiteSig', 'whiteMat', 'dewhiteMat'.\n");
	    }

	  /*Calculate the whitening */
	  whitenv (mixedsig, *(newE), *(newD), verbose, *(whitesig),
		   *(whiteningMatrix), *(dewhiteningMatrix));
	}
    }
  if (only > 2)
    {
      /*Calculating the ICA */
      
      /*Check some parameters
	The dimension of the data may have been reduced during PCA calculations.
	The original dimension is calculated from the data by default, and the
	number of IC is by default set to equal that dimension. */

      /*The number of IC's must be less or equal to the dimension of data */
      
      if (numOfIC > newDim)
	numOfIC = newDim;
      /*Show warning only if verbose = 'on' and user supplied a value for 'numOfIC' */
      if (verbose && (userNumOfIC == 1))
	{
	  fprintf(stderr, "Warning: estimating only %d independent components\n", numOfIC);
	  fprintf(stderr, "(Can't estimate more independent components than dimension of data)\n");
	}

      *(W) = mat_new_zeros (numOfIC, Dim);
      *(A) = mat_new_zeros (Dim, numOfIC);
      /*Calculate the ICA with fixed point algorithm. */
      
      *icasig = mat_new_zeros (numOfIC, NumOfSampl);
      
      fpICA (*(whitesig), *(whiteningMatrix), *(dewhiteningMatrix), approach,
	     numOfIC, g, finetune, a1, a2, myy, stabilization, epsilon,
	     maxNumIterations, maxFinetune, initState, guess, sampleSize,
	     verbose, *A, *W);
      
      /*Check for valid return */
      if (W)
	/*Add the mean back in. */
	{
	  if (verbose)
	    fprintf(stderr, "Adding the mean back to the data.\n");
	  
	  /*icasig = W * mixedsig + (W * mixedmean) * ones(1, NumOfSampl); */
	  mat_mul (*icasig, *W, mixedsig);
	  mat_temp_1 = mat_new (numOfIC, NumOfSampl);
	  vec_temp_1 = mat_vec_mul (*W, mixedmean);
	  
	  for (i = 0; i < NumOfSampl; i++)
	    mat_set_col (mat_temp_1, i, vec_temp_1);
	  vec_delete (vec_temp_1);
	  mat_add (*icasig, mat_temp_1);
	  mat_delete (mat_temp_1);
	  
	  /***icasig = W * mixedsig;***/
	}
      else
	mat_void (*icasig);	/*icasig = []; */
      
    }

  mat_delete (D);
  mat_delete (E);
  
  mat_delete (mixedsig);
  vec_delete (mixedmean);

  return 0;
}


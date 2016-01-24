/* Mex function for calculating cross-correlation of periodic signals
 *
 * Copyright (c) 2013, Jonathan Suever
 * Last Modified: 02-15-2013
 * Modified By: Suever (suever@gmail.com)
 */
#include "math.h"
#include "mex.h"
#include "blas.h"

#ifdef _MSC_VER
#include "inttypes.h"
#include "stdint.h"
#endif

/* Compute the standard deviation */
double compute_average(double *A, size_t nPoints)
{
    size_t i;
    double avg = 0.0;

    for (i=0; i<nPoints; i++){
        avg += A[i];
    }

    return avg / (double)nPoints;
}

/* Compute the standard deviation */
double compute_stdev(double *A, size_t nPoints, double avg)
{
    double diff;
    size_t i;
    double total=0.0;

    for (i=0; i<nPoints; i++){
        diff = A[i]-avg;
        total += diff * diff;
    }

    total /= (double)nPoints - 1.0;

    return sqrt(total);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A, *B, *C, *a, *b, *BMat, *inputs, avg_A, avg_B, std_A, std_B;
    double alpha, one = 1.0, zero = 0.0;

    int maxlag, minlag;
    INT32_T *lags, lag;

    size_t i, j, row, ind, nPoints, nLags, inc = 1;

    char *chn = "N";

    /* Input validation*/
    if (nrhs < 2)
        mexErrMsgTxt("You must provide two vectors as inputs");
    if (mxGetNumberOfElements(prhs[0]) != mxGetNumberOfElements(prhs[1]))
        mexErrMsgTxt("Input vectors must be the same length");

    /* Figure out the size of the inputs */
    nPoints = mxGetNumberOfElements(prhs[0]);

    /* Define lag range based upon the inputs */
    if (2 == nrhs) {                                    /* If no value supplied */
        maxlag = (nPoints + 1) / 2;
        minlag = -maxlag;
    } else {
        if (mxIsScalar(prhs[2])) {                      /* If scalar supplied */
            maxlag = (int)mxGetScalar(prhs[2]);
            minlag = -maxlag;
        } else {
            if (2 == mxGetNumberOfElements(prhs[2])) {  /* If vector supplied */
                inputs = mxGetPr(prhs[2]);
                maxlag = inputs[1];
                minlag = inputs[0];
            } else {
                mexErrMsgTxt("Third input should be either a scalar or two element array");
            }
        }
    }

    /* Make sure that for some reason we didn't end up with wonky lag times */
    if (maxlag <= minlag)
        mexErrMsgTxt("Your minimum lag should be less than your maximum lag");

    /* Grab the input data */
    A = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);

    a = (double*)mxMalloc(nPoints * sizeof(double));
    b = (double*)mxMalloc(nPoints * sizeof(double));

    /* Compute the mean of each of these guys */
    avg_A = compute_average(A, nPoints);
    avg_B = compute_average(B, nPoints);
    std_A = compute_stdev(A, nPoints, avg_A);
    std_B = compute_stdev(B, nPoints, avg_B);

    /* De-mean */
    for(i=0; i<nPoints; i++)
    {
        a[i] = A[i] - avg_A;
        b[i] = B[i] - avg_B;
    }

    /* Normalize by the standard deviation */
    alpha = 1.0/std_A;
    dscal(&nPoints, &alpha, a, &inc);

    alpha = 1.0 / std_B;
    dscal(&nPoints, &alpha, b, &inc);

    /* Compute the maximum lag (ceil) */
    nLags = maxlag - minlag + 1;

    /* Go ahead and create lags matrix */

    /* Go ahead and create the "B" matrix */
    BMat = (double*)mxMalloc(nPoints * nLags * sizeof(double));

    /* Create output variables */
    plhs[0] = mxCreateDoubleMatrix(nLags, 1, mxREAL);
    plhs[1] = mxCreateNumericMatrix(nLags, 1, mxINT32_CLASS, mxREAL);
    C       = mxGetPr(plhs[0]);
    lags    = (INT32_T*)mxGetPr(plhs[1]);

    row = 0;
    for(lag=minlag; lag<=maxlag; lag++)
    {
        for (i=0; i<nPoints; i++)
        {
            /* We have to add nPoints because we want to ensure % nPoints > 1
             * */
            ind = (i - lag + 2*nPoints) % nPoints;
            BMat[(i*nLags) + row] = b[ind];
        }
        /* Assign lags output variable value */
        lags[row] = lag;
        row++;
    }

    /* Perform matrix multiplication and normalization by (nPoints - 1) */
    alpha = 1.0 / (nPoints - 1.0);
    dgemm(chn, chn, &nLags, &inc, &nPoints, &alpha, BMat, &nLags, a, &nPoints, &zero, C, &nLags);

    /* Free up allocated values*/
    mxFree(BMat);
    mxFree(a);
    mxFree(b);
}

/* QT Clustering algorithm
 *
 *  Utilizes Euclidean distance to determine the largest possible clusters that
 *  meet the specified threshold value
 *
 * Copyright (c) 2012, Jonathan Suever
 * Last Modified: 06-05-2012
 * Modified By: Suever (suever@gmail.com)
 *
*/

#include "matrix.h"
#include "mex.h"
#include "math.h"
#include <stdio.h>
#include <string.h>
#include <float.h>

typedef struct {
    int member_a;
    int member_c;
    int clusterID;
} observation;

/* Compute the distance between all of the points initially */
void computeInitialDistances(double *G,double *D,mwSize nPoints,mwSize nDims) {
    /* Looping Variables */
    size_t i,j,d,col;

    double cumsum, diff, dist;

    for (i=0; i<nPoints; i++) {
        for (j=i+1; j<nPoints; j++) {
            cumsum = 0.0;
            for (d=0; d<nDims; d++) {
                col = d * nPoints;
                diff = G[i+col] - G[j+col];
                cumsum += diff * diff;
            }

            /* Mirror distances across the diagonal because the distance
             * between A and B is the same as the distance between B and A */
            dist = sqrt(cumsum);
            D[i + (j*nPoints)] = dist;
            D[j + (i*nPoints)] = dist;
        }
    }
}

/* Compute the RMSE between all of the points intially */
void computeInitialRMSE(double *G, double *D, mwSize nPoints, mwSize nDims) {
    /* Looping variables */
    size_t i,j,d,col;

    double cumsum, diff, rmse;

    for (i=0; i<nPoints; i++) {
        for (j=i+1; j<nPoints; j++) {
            cumsum = 0.0;
            for (d=0; d<nDims; d++) {
                col = d * nPoints;
                diff = G[i+col] - G[j+col];
                cumsum += diff * diff;
                /*tmp += pow(G[i+col] - G[j+col],2);*/
            }

            /* Normalize distance by the number of dataPoints */
            cumsum = cumsum / nDims;

            /* Mirror distances across the diagonal because the distance
             * between A and B is the same as the distance between B and A */
            rmse = sqrt(cumsum);
            D[i + (j*nPoints)] = rmse;
            D[j + (i*nPoints)] = rmse;
        }
    }
}

/* Actually perform the quality threshold clustering */
/*void qtCluster(double *D, double *qt, int *skip,
                            const size_t nPoints, double cutoff, int clustNdx) {*/
void qtCluster(double *D, observation *observations, const size_t nPoints, double cutoff, int clustNdx) {

    /* Booleans */
    int flag;

    /* Counters */
    size_t numRemaining, cardinality_a, offset, cardinality_c = 0;

    /* Diameter measurements */
    double diameter_a, diameter_c, minjdiam, tmp;

    observation *obj;

    double *jdiam;

    /* Loop variables */
    size_t i, j, k, closest_point_ndx;

    jdiam = (double*)mxMalloc(nPoints * sizeof(double));

    for (i=0; i<nPoints; i++) {
        /* Skip this point if it is already assigned */
        if (observations[i].clusterID)
            continue;

        flag = 1;
        cardinality_a = 1;
        diameter_a = 0.0;

        /* reset all member values to zeros */
        for (j=0; j<nPoints; j++)
            observations[j].member_a = 0;

        /* At least assign the current point to it's own cluster */
        observations[i].member_a = 1;

        while (flag & (cardinality_a < nPoints)) {
            for (j=0; j<nPoints; j++) {
                jdiam[j] = -1;

                obj = &observations[j];

                /* Ignore this point if we need to */
                if (obj->clusterID || obj->member_a)
                    continue;

                /* Precompute this for speed */
                offset = j * nPoints;

                /* Figure out the maximum distance from all member points */
                for (k=0; k<nPoints; k++) {
                    /* Skip if it isn't a number */
                    if ( observations[k].member_a == 0 || k == j)
                        continue;

                    tmp = D[k + offset];

                    /* Replace until we find the largest diameter */
                    if (tmp > jdiam[j])
                        jdiam[j] = tmp;
                }
            }

            /* After we have looped through all of the points we know the
             * maximum distance of each point from the cluster */

            minjdiam = DBL_MAX;

            /* Now find the smallest diameter */
            for (j=0; j<nPoints; j++) {
                if ((jdiam[j] < minjdiam) & (jdiam[j] > 0)) {
                    closest_point_ndx = j;
                    minjdiam = jdiam[j];
                }
            }

            /* Once we exceed the threshold, break out */
            if ((minjdiam > cutoff) | (diameter_a > cutoff))
                flag = 0;
            else {
                /* Keep adding successive points as long as we don't exceed */
                observations[closest_point_ndx].member_a = 1;
                cardinality_a += 1;

                /* Update cluster diameter if necessary */
                if (minjdiam > diameter_a)
                    diameter_a = minjdiam;
            }
        }

        /* If the cardinality of the current cluster is the largest */
        if (cardinality_a > cardinality_c) {
            for (j=0; j<nPoints; j++)
                observations[j].member_c = observations[j].member_a;
            cardinality_c = cardinality_a;
            diameter_c = diameter_a;
        }
    }

    numRemaining = 0;

    /* Remove added items from future consideration */
    for (j=0; j<nPoints; j++) {
        obj = &observations[j];
        if (obj->member_c){
            obj->clusterID = clustNdx;
        } else {
            /* Figure out how many points we have left to process */
            if (obj->clusterID == 0)
                numRemaining++;
        }
    }

    mxFree(jdiam);

    /* If we have to process more, then recurse to take care of them */
    if (numRemaining)
        qtCluster(D,observations,nPoints,cutoff,clustNdx+1);
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    /* Input variables */
    double *G, *D, *qtout, cutoff;

    int typeFlag = 0;   /* Specifies the initial. 0 = distance, 1 = RMSE */

    size_t nPoints, nDims;

    /* Loop Variables */
    size_t i;

    mwSize dims[1];

    observation *observations;
    observation default_observation = {0, 0, 0};

    /* Input checking */
    if (nrhs < 2) { mexErrMsgTxt("Must supply at least two inputs"); }

    if (nrhs == 3) typeFlag = (int)mxGetScalar(prhs[2]);

    if (mxGetM(prhs[0]) < 2) { mexErrMsgTxt("Data length must be > 1"); }
    if (nlhs > 1) { mexErrMsgTxt("Too many output arguments."); }

    /* Link up the inputs */
    G       = mxGetPr(prhs[0]);
    cutoff  = mxGetScalar(prhs[1]);
    nPoints = mxGetM(prhs[0]);
    nDims   = mxGetN(prhs[0]);
    dims[0] = nPoints;
    /* Link up outputs */
    plhs[0] = mxCreateNumericArray(1,dims,mxDOUBLE_CLASS,mxREAL);
    qtout = mxGetPr(plhs[0]);

    D = (double*)mxMalloc(nPoints*nPoints*sizeof(double));
    observations = (observation*)mxMalloc(nPoints*sizeof(observation));

    /* Initialize the struct */
    for (i=0; i<nPoints; i++) {
        observations[i] = default_observation;
    }

    /* Compute the initial distance matrix */
    if (typeFlag == 0) {
        computeInitialDistances(G,D,nPoints,nDims);
    } else if (typeFlag == 1) {
        computeInitialRMSE(G,D,nPoints,nDims);
    } else {
        mexErrMsgTxt("Invalid Type. Must be either 0 or 1");
    }

    /* Actually perform qt clustering */
    qtCluster(D, observations, nPoints, cutoff, 1);

    /* Save data into matlab output */
    for (i=0; i<nPoints; i++){
        qtout[i] = observations[i].clusterID;
    }

    mxFree(D);
    mxFree(observations);
}

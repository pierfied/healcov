//
// Created by pierfied on 10/19/19.
//

#ifndef HEALCOV_HEALCOV_H
#define HEALCOV_HEALCOV_H

typedef struct {
    long nside;
    long npix;
    long *mask_inds;
    long tree_depth;
    long nsamps;
    double *theta_samps;
    double *xi_samps;
} HealcovArgs;

double *build_cov(HealcovArgs args);

#endif //HEALCOV_HEALCOV_H

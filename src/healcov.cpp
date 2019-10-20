#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
//
// Created by pierfied on 10/19/19.
//

#include "healcov.h"
#include <healpix_cxx/healpix_base.h>
#include <cmath>

double *build_cov(HealcovArgs args) {
    long nside = args.nside;
    long npix = args.npix;
    long *mask_inds = args.mask_inds;
    long tree_depth = args.tree_depth;
    long sub_npix = tree_depth * tree_depth;
    long nsamps = args.nsamps;
    double *theta_samps = args.theta_samps;
    double *xi_samps = args.xi_samps;

    long sub_nside = tree_depth * nside;
    const nside_dummy dummy;
    T_Healpix_Base<long> sub_base = T_Healpix_Base<long>(sub_nside, NEST, dummy);

    double sub_vecs_x[npix * sub_npix];
    double sub_vecs_y[npix * sub_npix];
    double sub_vecs_z[npix * sub_npix];
#pragma omp parallel for
    for (long i = 0; i < npix; ++i) {
        long sub_ind0 = i * sub_npix;
        long sub_pix0 = mask_inds[i] * sub_npix;

        for (long j = 0; j < sub_npix; ++j) {
            vec3 vec_ij = sub_base.pix2vec(sub_pix0 + j);

            long sub_ind = sub_ind0 + j;

            sub_vecs_x[sub_ind] = vec_ij.x;
            sub_vecs_y[sub_ind] = vec_ij.y;
            sub_vecs_z[sub_ind] = vec_ij.z;
        }
    }

    double theta0 = theta_samps[0];
    double dtheta = theta_samps[1] - theta_samps[0];
    double *cov = new double[npix * npix];
#pragma omp parallel for
    for (long i = 0; i < npix; ++i) {
        long sub_ind_i = i * sub_npix;

        for (long j = 0; j < npix; ++j) {
            long sub_ind_j = j * sub_npix;

            double cov_ij = 0;
            for (long k = sub_ind_i; k < sub_ind_i + sub_npix; ++k) {
                for (long l = sub_ind_j; l < sub_ind_j + sub_npix; ++l) {
                    double sep = acos(std::max<double>(std::min<double>(
                            sub_vecs_x[k] * sub_vecs_x[l] + sub_vecs_y[k] * sub_vecs_y[l] +
                            sub_vecs_z[k] * sub_vecs_z[l], 1), -1));

                    long theta_ind = (sep - theta0) / dtheta;
                    double frac_low = (sep - theta_samps[theta_ind]) / dtheta;

                    cov_ij += frac_low * xi_samps[theta_ind] + (1 - frac_low) * xi_samps[theta_ind + 1];
                }
            }

            cov[i * npix + j] = cov_ij / (sub_npix * sub_npix);
        }
    }

    return cov;
}

#pragma clang diagnostic pop
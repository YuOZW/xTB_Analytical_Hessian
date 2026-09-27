/* This file is part of tblite.
 * SPDX-Identifier: LGPL-3.0-or-later
 */

#include <math.h>
#include <stdio.h>
#include "tblite.h"

static int print_error(tblite_error error, tblite_context ctx)
{
    char message[512] = {0};
    if (error && tblite_check(error)) {
        tblite_get_error(error, message, NULL);
        fprintf(stderr, "[Fatal] %s\n", message);
        return 1;
    }
    if (ctx && tblite_check(ctx)) {
        tblite_get_context_error(ctx, message, NULL);
        fprintf(stderr, "[Fatal] %s\n", message);
        return 1;
    }
    return 0;
}

int main(void)
{
    tblite_error error = tblite_new_error();
    tblite_context ctx = tblite_new_context();
    tblite_result res = tblite_new_result();
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    const int nat = 3;
    int num[3] = {8, 1, 1};
    double xyz[9] = {
         0.000000000000, 0.000000000000, 0.0,
         1.432336730000, 1.107152660000, 0.0,
        -1.432336730000, 1.107152660000, 0.0
    };
    double charge = 0.0;
    int uhf = 0;
    double hessian[81] = {0.0};
    double frequencies[3] = {0.0};
    int nvib = 0, nrigid = 0, nimag = 0;
    double max_asym = 0.0, max_trans = 0.0;
    int i, j, a, c;
    int status = 1;

    tblite_set_context_verbosity(ctx, 0);
    mol = tblite_new_structure(error, nat, num, xyz, &charge, &uhf, NULL, NULL);
    if (print_error(error, ctx)) goto cleanup;
    calc = tblite_new_gfn1_calculator(ctx, mol, NULL);
    if (!calc || print_error(error, ctx)) goto cleanup;
    tblite_set_calculator_temperature(ctx, calc, 0.0);
    if (print_error(error, ctx)) goto cleanup;

    tblite_get_hessian(ctx, mol, calc, res);
    if (print_error(error, ctx)) goto cleanup;
    tblite_get_result_hessian(error, res, hessian);
    tblite_get_result_number_of_vibrations(error, res, &nvib);
    tblite_get_result_number_of_rigid_modes(error, res, &nrigid);
    tblite_get_result_number_of_imaginary_modes(error, res, &nimag);
    tblite_get_result_frequencies(error, res, frequencies);
    if (print_error(error, ctx)) goto cleanup;

    if (nvib != 3 || nrigid != 6 || nimag < 0) goto cleanup;
    for (i = 0; i < 3; ++i) if (!isfinite(frequencies[i])) goto cleanup;
    for (i = 0; i < 9; ++i) {
        for (j = 0; j < 9; ++j) {
            double d;
            if (!isfinite(hessian[i*9+j])) goto cleanup;
            d = fabs(hessian[i*9+j] - hessian[j*9+i]);
            if (d > max_asym) max_asym = d;
        }
    }
    if (max_asym > 1.0e-10) goto cleanup;

    for (i = 0; i < 9; ++i) {
        for (c = 0; c < 3; ++c) {
            double sum = 0.0;
            for (a = 0; a < nat; ++a) sum += hessian[i*9 + 3*a + c];
            if (fabs(sum) > max_trans) max_trans = fabs(sum);
        }
    }
    if (max_trans > 1.0e-8) goto cleanup;

    status = 0;

cleanup:
    if (status != 0) print_error(error, ctx);
    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res);
    return status;
}

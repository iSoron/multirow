/*
 * Copyright (C) 2016 Álinson Santos Xavier <isoron@gmail.com>
 *
 * This file is part of Loop Habit Tracker.
 *
 * Loop Habit Tracker is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * Loop Habit Tracker is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include<cblas.h>
#include<lapacke.h>
#include <multirow/util.h>

void LINALG_scale(int n, double *x, double alpha)
{
    cblas_dscal(n, alpha, x, 1);
}

double LINALG_dot(int n, double *x, double *y)
{
    return cblas_ddot(n, x, 1, y, 1);
}

double LINALG_norm(int n, double *x)
{
    return cblas_dasum(n, x, 1);
}

int LINALG_solve(int n, double *A, double *b, double *x)
{
    int rval = 0;
    int ignored[n];
    double A_copy[n * n];

    memcpy(x, b, n * sizeof(double));
    memcpy(A_copy, A, n * n * sizeof(double));

    rval = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, 1, A_copy, n, ignored, x, 1);
    abort_if(rval, "LAPACKE_dgesv failed");

CLEANUP:
    return rval;
}
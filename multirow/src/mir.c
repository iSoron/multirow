/* Copyright (c) 2015 Alinson Xavier
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <multirow/lp.h>
#include <multirow/mir.h>
#include <multirow/util.h>

static double h(double a)
{
    return fmax(a, 0);
}

static double f(double a,
                double b)
{
    return frac(b) * floor(a) + fmin(frac(a), frac(b));
}

int MIR_generate_cut(const struct Row *row,
                     const char *column_types,
                     struct Row *cut)
{
    int rval = 0;

    cut->pi = (double*) malloc(row->nz * sizeof(double));
    cut->indices = (int*) malloc(row->nz * sizeof(int));
    abort_if(!cut->pi, "could not allocate cut->pi");
    abort_if(!cut->indices, "could not allocate cut->indices");

    cut->head = 0;
    cut->nz = row->nz;
    cut->pi_zero = -frac(row->pi_zero) * ceil(row->pi_zero);

    for (int i = 0; i < cut->nz; i++)
    {
        int idx = row->indices[i];
        double value = row->pi[i];

        cut->indices[i] = idx;

        if (column_types[idx] == MILP_INTEGER)
            cut->pi[i] = -f(value, row->pi_zero);
        else
            cut->pi[i] = -h(value);
    }

CLEANUP:
    return rval;
}

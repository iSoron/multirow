/* Copyright (c) 2015-2017 Alinson Xavier
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

#pragma once

/**
 * Receives an n-dimensional vector x, a scalar alpha and sets x <- alpha * x
 */
void LINALG_scale(int n, double *x, double alpha);

/**
 * Receives two n-dimensional vectors x and y and returns the dot product of x
 * and y.
 */
double LINALG_dot(int n, double *x, double *y);

/**
 * Receives an n-dimensional vector x and returns the 1-norm of x.
 */
double LINALG_norm(int n, double *x);

/**
 * Given a full rank m-by-n matrix A and an m-dimensional vector b, this
 * function finds x such that Ax = b. Returns zero if the operation is
 * successful and non-zero otherwise.
 */
int LINALG_solve(int n, int m, double *A, double *b, double *x);
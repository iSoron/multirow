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
#ifndef MULTIROW_GREEDY_ND_H
#define MULTIROW_GREEDY_ND_H

int GREEDY_ND_next_lattice_point(int dim,
                                 const double *lb,
                                 const double *ub,
                                 double *p,
                                 int *finished);

int GREEDY_create_psi_lp(const int nrows,
                         const int nrays,
                         const double *f,
                         const double *rays,
                         const double *beta,
                         struct LP *lp);

int GREEDY_ND_psi(const int nrows,
                  const int nrays,
                  const double *f,
                  const double *rays,
                  const double *beta,
                  const double *q,
                  const double q_scale,
                  struct LP *lp,
                  double *value);

int GREEDY_ND_pi(const int nrows,
                  const int nrays,
                  const double *f,
                  const double *rays,
                  const double *beta,
                  const double *q,
                  const double q_scale,
                  struct LP *lp,
                  double *value);

int GREEDY_ND_generate_cut(int nrows,
                           int nrays,
                           const double *f,
                           const double *rays,
                           double *beta);

int GREEDY_ND_bound(int nrows,
                    int nrays,
                    const double *f,
                    const double *rays,
                    const double *x,
                    const double *beta,
                    double *epsilon,
                    int *tx);

int GREEDY_ND_cone_bound(int nrows,
                         int nrays,
                         const double *f,
                         const double *rays,
                         const int *rx,
                         const double *x,
                         const double *beta,
                         double *epsilon);

int GREEDY_ND_find_violated_cone(int nrows,
                                 int nrays,
                                 const double *f,
                                 const double *rays,
                                 const double *x,
                                 const double *beta,
                                 double epsilon,
                                 int *rx,
                                 double *sbar,
                                 int *violated_found);

int GREEDY_ND_find_tight_rays(int nrows,
                              int nrays,
                              const double *f,
                              const double *rays,
                              const double *x,
                              const double *beta,
                              double epsilon,
                              int *tx);

int GREEDY_ND_scale_to_ahull(int nrows,
                             int nrays,
                             const double *rays,
                             const int *rx,
                             const double *beta,
                             double epsilon,
                             const double *d,
                             double *alpha);

#endif //MULTIROW_GREEDY_ND_H

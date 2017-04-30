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
#ifndef MULTIROW_INFINITY_ND_H
#define MULTIROW_INFINITY_ND_H

int INFINITY_create_psi_lp(const struct ConvLFreeSet *lfree, struct LP *lp);

int INFINITY_psi(const int nrows,
                 const double *q,
                 const double q_scale,
                 struct LP *lp,
                 double *value);

int INFINITY_pi(const int nrows,
                const double *q,
                const double q_scale,
                struct LP *lp,
                double *value);

int INFINITY_ND_generate_lfree(const struct MultiRowModel *model,
                               struct ConvLFreeSet *lfree);

#endif //MULTIROW_INFINITY_ND_H

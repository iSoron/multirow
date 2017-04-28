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
#ifndef MULTIROW_GREEDY_H
#define MULTIROW_GREEDY_H

int GREEDY_write_sage_file(int nrows,
                           int nrays,
                           const double *f,
                           const double *rays,
                           const double *bounds,
                           const char *filename);

int GREEDY_generate_cut(int nrows,
                        struct Row **rows,
                        const char *column_types,
                        struct Row *cut);

#endif //MULTIROW_GREEDY_H

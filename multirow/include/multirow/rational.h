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

#ifndef MULTIROW_RATIONAL_H
#define MULTIROW_RATIONAL_H

struct _rational
{
    long num;
    unsigned long den;
};

typedef struct _rational Rational[1];

long gcd(long a,
         long b);

long lcm(long a,
         long b);

#endif //MULTIROW_RATIONAL_H

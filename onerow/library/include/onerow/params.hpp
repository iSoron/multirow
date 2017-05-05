/*
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

#ifndef PARAMS_HPP_
#define PARAMS_HPP_

const double ZERO_CUTOFF = 1e-8;

const double MAX_CUT_DYNAMISM  = 1e6;
const double MIN_CUT_VIOLATION = 1e-6;

const long REDUCE_FACTOR_RHS          = 1000000;
const long REDUCE_FACTOR_R1           = 1000;
const long REDUCE_FACTOR_COEFFICIENT  = 1000000;

const int MAX_R1_RAYS   = 1000000;
const int MAX_CUT_DEPTH = 1000000;
const int MAX_GOOD_ROWS = 1000000;

const int ETA_UPDATE_INTERVAL = 300;
const unsigned int MAX_CUT_BUFFER_SIZE = 100;

#define INTERSECTION_CUT_USE_DOUBLE
// #define ENABLE_EXTENDED_STATISTICS
// #define PRETEND_TO_ADD_CUTS

#endif /* PARAMS_HPP_ */

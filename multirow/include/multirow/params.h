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

#ifndef PROJECT_PARAMS_H
#define PROJECT_PARAMS_H

/*
 * Error margin for floating point comparisons.
 */
#define EPSILON 1e-8

/*
 * Available log levels, in decreasing level of verboseness, are:
 *
 *      LOG_LEVEL_VERBOSE
 *      LOG_LEVEL_DEBUG
 *      LOG_LEVEL_INFO
 *      LOG_LEVEL_WARNING
 *      LOG_LEVEL_ERROR
 */
#define LOG_LEVEL LOG_LEVEL_INFO

/*
 * Maximum bounding-box size for naive algorithm
 */
#define MAX_BOX_SIZE 10000

/*
 * Maximum number of sets that should be considered
 */
#define MAX_N_SETS 1000

/*
 * Number of rays that should be generated per set.
 */
#define N_RAYS 100

#define ONLY_CUT -1

/*
 * Time limit for the computation (user time, in seconds).
 */
#define MAX_TOTAL_TIME -1

/*
 * Time limit for each CPLEX LP or MIP (in seconds).
 */
#define CPLEX_TIMEOUT 1.0
#define CG_TIMEOUT 900

#define MAX_N_RAYS 100
#define MAX_SELECTED_COMBINATIONS 5000
#define MAX_SELECTED_ROWS 300
#define MAX_LATTICE_POINTS 100000

#define INFINITY_BIG_E 1024

#define MAX_CUT_DYNAMISM 1e8
#define INTEGRALITY_THRESHOLD 0.49

extern int BOOST_VAR;
extern double BOOST_FACTOR;
extern int SHOULD_DUMP_CUTS;
extern int DUMP_CUT_N;

extern int ENABLE_LIFTING;
extern int MIN_N_ROWS;
extern int MAX_N_ROWS;

#define ERR_NO_CUT 2
#define ERR_MIP_TIMEOUT 3

#endif //PROJECT_PARAMS_H

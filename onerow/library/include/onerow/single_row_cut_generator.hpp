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

#ifndef SINGLE_ROW_CUT_GENERATOR_HPP_
#define SINGLE_ROW_CUT_GENERATOR_HPP_

#include <ilcplex/cplex.h>
#include <vector>
#include "geometry.hpp"
using std::vector;

/**
 * Models a linear constraint.
 */
struct Constraint {

	/**
	 * Vector that holds the coefficients of the variables.
	 */
	svec pi;

	/**
	 * The right hand side of the constraint.
	 */
	rational pi_zero;

	int depth;

	/**
	 * Comparator used to sort the constraints.
	 *
	 * @returns bool True if this constraint should come before the given
	 * constraint when sorting.
	 */
	bool operator<(const Constraint &other) const;

	bool operator==(const Constraint &other) const;

	rational get_violation(const rational *x);
};

/**
 * Models a single row from the simplex tableau.
 */
struct Row {
	/**
	 * Constraint
	 */
	Constraint c;

	/**
	 * Index of the basic variable this row corresponds to.
	 */
	int basic_var_index;

	bool* is_integer;
	double* reduced_costs;
	double cost_cutoff;
};

/**
 * A single row cut generator receives a row from the simplex tableau and generates one
 * or more cuts that invalidate the current basic solution.
 */
class SingleRowCutGenerator {
protected:
	const Row& row;

public:

	/**
	 * Constructs a new generator that will generate cuts from the provided tableau row.
	 */
	SingleRowCutGenerator(Row &r);

	/**
	 * Destructor.
	 */
	virtual ~SingleRowCutGenerator() {};

	/**
	 * Returns true if more cuts can be generated from the tableau row.
	 */
	virtual bool has_next() = 0;

	/**
	 * Retrieves a cut generated from the tableau row.
	 *
	 * @throws std::out_of_bounds if no more cuts can be generated from the current tableau row.
	 */
	virtual Constraint* next() = 0;
};

std::ostream& operator<<(std::ostream& os, const Constraint &c);
std::ostream& operator<<(std::ostream& os, const Row &r);

#endif

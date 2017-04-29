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

#ifndef MIR_CUT_GENERATOR_HPP_
#define MIR_CUT_GENERATOR_HPP_

#include "single_row_cut_generator.hpp"

class MIRCutGenerator: public SingleRowCutGenerator {
private:
	bool finished;
	rational f(rational q1, rational q2);
	rational h(rational q);

public:
	MIRCutGenerator(Row &row);
	~MIRCutGenerator();

	bool has_next();
	Constraint* next();
};

#endif /* MIR_CUT_GENERATOR_HPP_ */

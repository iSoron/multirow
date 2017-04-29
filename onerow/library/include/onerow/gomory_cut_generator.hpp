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

#ifndef GOMORY_CUT_GENERATOR_HPP_
#define GOMORY_CUT_GENERATOR_HPP_

#include "single_row_cut_generator.hpp"

/**
 * This class can be used to generate classic (not fractional) Gomory cuts.
 * The cuts are only valid for models with integral, non-negative variables.
 */
class GomoryCutGenerator: public SingleRowCutGenerator {
private:
	bool finished;

public:
	GomoryCutGenerator(Row &row);
	~GomoryCutGenerator();

	bool has_next();
	Constraint* next();
};

#endif /* GOMORY_CUT_GENERATOR_HPP_ */

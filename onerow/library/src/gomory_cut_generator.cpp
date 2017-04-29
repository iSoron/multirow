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

#include <stdexcept>
#include <onerow/gomory_cut_generator.hpp>

GomoryCutGenerator::GomoryCutGenerator(Row &r) :
		SingleRowCutGenerator(r), finished(false)
{

}

GomoryCutGenerator::~GomoryCutGenerator()
{

}

bool GomoryCutGenerator::has_next()
{
	return !finished;
}

Constraint* GomoryCutGenerator::next()
{
	if (!has_next())
		throw std::out_of_range("");

	Constraint *cut = new Constraint;
	int nz = row.c.pi.nz();
	cut->pi_zero = row.c.pi_zero.floor();

	cut->pi.resize(row.c.pi.size());
	for (int i = 0; i < nz; i++)
		cut->pi.push(row.c.pi.index(i), row.c.pi.value(i).floor());

	finished = true;
	return cut;
}

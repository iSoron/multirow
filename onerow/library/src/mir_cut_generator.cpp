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
#include <onerow/stats.hpp>
#include <onerow/mir_cut_generator.hpp>

MIRCutGenerator::MIRCutGenerator(Row &r) :
		SingleRowCutGenerator(r), finished(false)
{

}

MIRCutGenerator::~MIRCutGenerator()
{

}

bool MIRCutGenerator::has_next()
{
	return !finished;
}

rational MIRCutGenerator::h(rational a)
{
	if (a > 0) return a;
	else return 0;
}

rational MIRCutGenerator::f(rational a, rational b)
{
	if (a.frac() <= b.frac())
		return b.frac() * a.floor() + a.frac();
	else
		return b.frac() * a.ceil();
}

Constraint* MIRCutGenerator::next()
{
	if (!has_next())
		throw std::out_of_range("");

	Constraint *cut = new Constraint;
	int nz = row.c.pi.nz();
	cut->pi.resize(row.c.pi.size());
	cut->pi_zero = -row.c.pi_zero.frac() * row.c.pi_zero.ceil();

	for (int i = 0; i < nz; i++)
	{
		int idx = row.c.pi.index(i);

		if (row.is_integer[idx])
			cut->pi.push(idx, -f(row.c.pi.value(i), row.c.pi_zero));
		else
			cut->pi.push(idx, -h(row.c.pi.value(i)));
	}

	cut->depth = 0;

	finished = true;
	return cut;
}

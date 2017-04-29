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

#include <onerow/single_row_cut_generator.hpp>

SingleRowCutGenerator::SingleRowCutGenerator(Row &r) : row(r)
{

}

#define return_if_neq(a,b) if((a)<(b)) return true; if((a)>(b)) return false;

bool Constraint::operator<(const Constraint &other) const
{
	return_if_neq(pi.nz(), other.pi.nz());
	return_if_neq(pi_zero, other.pi_zero);

	int nz = pi.nz();
	for (int i = 0; i < nz; i++)
	{
		return_if_neq(pi.index(i), other.pi.index(i));
		return_if_neq(pi.value(i), other.pi.value(i));
	}

	return false;
}

bool Constraint::operator==(const Constraint &other) const
{
	return !(operator<(other) || other.operator<(*this));
}

rational Constraint::get_violation(const rational *x)
{
	rational v(0);

	int nz = pi.nz();
	for (int i = 0; i < nz; i++)
		v += pi.value(i) * x[pi.index(i)];
	v -= pi_zero;

	return v;
}

std::ostream& operator<<(std::ostream& os, const Constraint &c)
{
	int nz = c.pi.nz();
	for (int k = 0; k < nz; k++)
	{
		os << c.pi.value(k) << " x" << c.pi.index(k) << " ";
	}
	os << "<= " << c.pi_zero;
	return os;
}

std::ostream& operator<<(std::ostream& os, const Row &r)
{
	os << "[" << r.basic_var_index << "] " << r.c;
	return os;
}

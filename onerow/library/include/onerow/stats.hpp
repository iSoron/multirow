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

#ifndef STATS_HPP_
#define STATS_HPP_

#include <string>
using std::string;

namespace Stats {
	void init();
	void add_cut(int depth);
	void add_generated_cut(int round, int depth);
	void set_solution(int round, double sol, string status);
	void set_input_filename(string n);

	void add_trivial_lifting_m(unsigned long m);
	void add_coefficient(bool integral);

	void write_stats(string filename);

	void start_timer();
	void end_timer();
};

double get_current_time(void);
void time_printf(const char *fmt, ...);

#endif /* STATS_HPP_ */

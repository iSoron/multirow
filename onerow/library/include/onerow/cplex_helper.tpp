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

#ifndef CPLEX_HELPER_TPP_
#define CPLEX_HELPER_TPP_
#include <cmath>
#include <omp.h>
#include <ctime>
#include <thread>
#include <unistd.h>

#include <onerow/stats.hpp>
#include <onerow/cplex_helper.hpp>
#include <onerow/params.hpp>

using std::cout;
using std::endl;


template<class Generator>
int CplexHelper::add_single_row_cuts(int max_rows)
{
	total_cuts = 0;

	if(n_good_rows > 0)
	{
		n_good_rows = 0;
		delete good_rows;
	}

	find_good_rows(max_rows);

	eta_reset();
	eta_count = 0;
	eta_total = n_good_rows;
	std::thread eta(&CplexHelper::eta_print, this);

	Stats::start_timer();

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < n_good_rows; i++)
	{
		Row *row = get_tableau_row(good_rows[i]);

		Generator generator(*row);

		while (generator.has_next())
		{
			Constraint *cut = generator.next();

			if (cut->pi.nz() == 0)
			{
				delete cut;
				continue;
			}

			#ifdef PRETEND_TO_ADD_CUTS

				delete(cut);

			#else

				#pragma omp critical(cplex)
				{
					add_cut(cut);
				}

			#endif

			#ifdef ENABLE_EXTENDED_STATISTICS
				Stats::add_generated_cut(current_round, cut->depth);
			#endif
		}

		#pragma omp critical
		{
			eta_count++;
		}

		delete row;
	}

	Stats::end_timer();

	eta.join();

	flush_cuts();
	time_printf("Added %d violated cuts...\n", total_cuts);

	return 0;
}

#endif /* CPLEX_HELPER_TPP_ */

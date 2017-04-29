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

#include <cstdio>
#include <cmath>
#include <cstdarg>
#include <cassert>
#include <sys/resource.h>
#include <sys/time.h>
#include <iostream>
#include <iomanip>

#include <onerow/params.hpp>
#include <onerow/stats.hpp>

namespace Stats
{
	const int MAX_ROUNDS = 100;
	const int MAX_TIMERS = 100;

	string input_filename;

	double opt_value[MAX_ROUNDS];
	string opt_sol_status[MAX_ROUNDS];

	unsigned long n_cuts_total = 0;
	unsigned long n_cuts_depth[MAX_CUT_DEPTH] = { 0 };

	unsigned long n_generated_cuts_total = 0;
	unsigned long n_generated_cuts_round[MAX_CUT_DEPTH] = { 0 };
	unsigned long n_generated_cuts_depth[MAX_CUT_DEPTH] = { 0 };

	unsigned long trivial_lifting_m_count = 0;
	unsigned long trivial_lifting_m_sum = 0;
	unsigned long trivial_lifting_m_max = 0;

	unsigned long n_coefficients = 0;
	unsigned long n_integral_coefficients = 0;

	int n_timers = 0;
	double current_timer_start;
	double timers[MAX_TIMERS] = {0};

	void init()
	{
		input_filename = "";

		for (int i = 0; i < MAX_ROUNDS; i++)
		{
			opt_value[i] = INFINITY;
			opt_sol_status[i] = "";
		}
	}

	void add_cut(int depth)
	{
		n_cuts_total++;
		n_cuts_depth[depth]++;
	}

	void add_generated_cut(int round, int depth)
	{
		n_generated_cuts_total++;
		n_generated_cuts_round[round]++;
		n_generated_cuts_depth[depth]++;
	}

	void add_trivial_lifting_m(unsigned long m)
	{
		trivial_lifting_m_count++;
		trivial_lifting_m_sum += m;
		trivial_lifting_m_max = std::max(trivial_lifting_m_max, m);
	}

	void add_coefficient(bool integral)
	{
		n_coefficients++;
		if(integral) n_integral_coefficients++;
	}

	void set_solution(int round, double sol, string status)
	{
		opt_value[round] = sol;
		opt_sol_status[round] = status;
	}

	void set_input_filename(string n)
	{
		input_filename = n;
	}

	void write_stats(string f)
	{
		FILE *out = fopen(f.c_str(), "w");

		fprintf(out, "input_file:\n  %s\n", input_filename.c_str());

		// solution value and status
		fprintf(out, "sol_value:\n");
		for (int i = 0; i < MAX_ROUNDS; i++)
			if (opt_value[i] < INFINITY)
				fprintf(out, "  %d: %-.6f\n", i, opt_value[i]);

		fprintf(out, "sol_status:\n");
		for (int i = 0; i < MAX_ROUNDS; i++)
			if (opt_value[i] < INFINITY)
				fprintf(out, "  %d: %s\n", i, opt_sol_status[i].c_str());

		// added cuts
		if(n_cuts_total > 0)
		{
			fprintf(out, "n_added_cuts:\n");

			fprintf(out, "  total:\n    %ld\n", n_cuts_total);

			fprintf(out, "  depth:\n");
			for (int i = 0; i < MAX_CUT_DEPTH; i++)
				if (n_cuts_depth[i] > 0)
					fprintf(out, "    %d: %ld\n", i, n_cuts_depth[i]);
		}

		// generated cuts
		if(n_generated_cuts_total > 0)
		{
			fprintf(out, "n_generated_cuts:\n");

			fprintf(out, "  total:\n    %ld\n", n_generated_cuts_total);

			fprintf(out, "  round:\n");
				for (int i = 0; i < MAX_ROUNDS; i++)
					if (n_generated_cuts_round[i] > 0)
						fprintf(out, "    %d: %ld\n", i, n_generated_cuts_round[i]);

			fprintf(out, "  depth:\n");
			for (int i = 0; i < MAX_CUT_DEPTH; i++)
				if (n_generated_cuts_depth[i] > 0)
					fprintf(out, "    %d: %ld\n", i, n_generated_cuts_depth[i]);
		}

		// trivial lifting
		if(trivial_lifting_m_count > 0)
		{
			double integral_coefficients = ((double) n_integral_coefficients) / n_coefficients;
			double average_m = ((double) (trivial_lifting_m_sum))
					/ trivial_lifting_m_count;
			double slowdown = integral_coefficients * average_m;

			fprintf(out, "trivial_lifting:\n");
			fprintf(out, "  max_m: %ld\n", trivial_lifting_m_max);
			fprintf(out, "  average_m: %.6lf\n", average_m);
			fprintf(out, "  integral_coefficients: %.6lf\n", integral_coefficients);
			fprintf(out, "  slowdown: %.6lf\n", slowdown);
		}

		if(n_timers > 0)
		{
			fprintf(out, "timers:\n");
			for(int i=0; i<n_timers; i++)
				fprintf(out, "  %d: %.4lf\n", i+1, timers[i]);
		}

		fprintf(out, "cut_speed:\n");
		for(int i=0; i<n_timers; i++)
			fprintf(out, "  round_%d: %.4lf\n", i+1, timers[i] / n_generated_cuts_round[i+1]);

	}


	void start_timer()
	{
		assert(current_timer_start == 0);

		time_printf("Starting timer %d...\n", n_timers+1);
		current_timer_start = get_current_time();
		timers[n_timers] = 0;
	}

	void end_timer()
	{
		timers[n_timers] = get_current_time() - current_timer_start;

		time_printf("Ending timer %d: %.2lfs\n", n_timers+1, timers[n_timers]);

		current_timer_start = 0;
		n_timers++;
	}


}

double get_current_time(void)
{
	struct rusage ru;

	getrusage (RUSAGE_SELF, &ru);

	return ((double) ru.ru_utime.tv_sec) +
		   ((double) ru.ru_utime.tv_usec)/1000001.0;
}

double initial_time = 0;

void time_printf(const char *fmt, ...)
{
	if(initial_time == 0)
		initial_time = get_current_time();

    std::cout << std::setw(40) << Stats::input_filename << " ";
	printf("[%10.2lf] ", get_current_time() - initial_time);

	va_list args;
	va_start(args, fmt);
	vprintf(fmt, args);
	va_end(args);

	fflush(stdout);
}

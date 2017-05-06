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
#include <cstdlib>
#include <string>
#include <unistd.h>
#include <vector>
#include <ilcplex/cplex.h>

#include <onerow/cplex_helper.hpp>
#include <onerow/cplex_helper.tpp>
#include <onerow/stats.hpp>

#include <onerow/gomory_cut_generator.hpp>
#include <onerow/mir_cut_generator.hpp>
#include <onerow/wedge_cut_generator.hpp>

using namespace std;

char *input_filename = 0;
char *stats_filename = 0;
char *sol_filename = 0;
char *basis_filename = 0;
bool enable_gomory_cuts = false;
bool enable_wedge_cuts = false;
bool enable_mir_cuts = false;

int c;
extern char *optarg;
char usage[] = "usage: %s [-gmw] -f model.mps [-b basis.bas] [-s stats.yaml]\n";

void read_params(int argc, char **argv)
{
	while ((c = getopt(argc, argv, "gwmf:s:x:b:")) != -1)
	{
		switch (c)
		{
		case 'g':
			enable_gomory_cuts = true;
			break;
		case 'w':
			enable_wedge_cuts = true;
			break;
		case 'm':
			enable_mir_cuts = true;
			break;
		case 'f':
			input_filename = optarg;
			break;
		case 'b':
			basis_filename = optarg;
			break;
		case 's':
			stats_filename = optarg;
			break;
        case 'x':
			sol_filename = optarg;
			break;
		}
	}

	if (!input_filename)
	{
		fprintf(stderr, "%s: missing model filename\n", argv[0]);
		fprintf(stderr, usage, argv[0]);
		exit(1);
	}
}

int main(int argc, char **argv)
{
	read_params(argc, argv);

	Stats::init();

	int status;
	CPXENVptr env = CPXopenCPLEX(&status);
	CPXLPptr lp = CPXcreateprob(env, &status, "");
	CplexHelper cplexHelper(env, lp);

	CPXsetlogfile(env, NULL);
	CPXsetintparam(env, CPX_PARAM_PREIND, CPX_OFF); // disable presolve
	CPXsetintparam(env, CPX_PARAM_DATACHECK, CPX_ON); // check consistency
	CPXsetintparam(env, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON); // numerical precision

	Stats::set_input_filename(string(input_filename));

	time_printf("Using OpenMP (%d threads)\n", omp_get_max_threads());

	// reads input file
	time_printf("Reading input file: %s...\n", input_filename);
	status = CPXreadcopyprob(env, lp, input_filename, NULL);
	if (status)
	{
		fprintf(stderr, "could not read input file (%d)\n", status);
		return 1;
	}

	cplexHelper.read_columns();

	// read solution
	if(sol_filename)
	{
		time_printf("Reading solution file: %s...\n", sol_filename);
		FILE *f = fopen(sol_filename, "r");
		if(!f)
		{
			fprintf(stderr, "Could not open solution file (%s).", sol_filename);
			return 1;
		}

		double *solution = new double[cplexHelper.n_cols];
		for(int i=0; i<cplexHelper.n_cols; i++)
			fscanf(f, "%lf", &solution[i]);

		cplexHelper.optimal_solution = solution;
	}


	// relaxes integrality
	CPXchgprobtype(env, lp, CPXPROB_LP);

	if (basis_filename)
	{
		time_printf("Loading basis from %s...\n", basis_filename);
		status = CPXreadcopybase(env, lp, basis_filename);
		if(status)
		{
			fprintf(stderr, "could not read basis file");
			return 1;
		}
	}

	time_printf("Solving first relaxation...\n");
	cplexHelper.solve(true);


	if (enable_gomory_cuts)
	{
		time_printf("Generating Gomory cuts...\n");
		cplexHelper.add_single_row_cuts<GomoryCutGenerator>(0);
		cplexHelper.solve(true);
	}

	if (enable_mir_cuts)
	{
		time_printf("Generating MIR cuts...\n");
		cplexHelper.add_single_row_cuts<MIRCutGenerator>(0);
		cplexHelper.solve(true);
	}

	if (enable_wedge_cuts)
	{
		time_printf("Generating wedge cuts...\n");
		cplexHelper.add_single_row_cuts<WedgeCutGenerator>(MAX_GOOD_ROWS);
		cplexHelper.solve(true);
	}

	if(stats_filename != 0)
	{
		time_printf("Writting stats: %s...\n", stats_filename);
		Stats::write_stats(string(stats_filename));
	}

	time_printf("Done.\n");
	CPXcloseCPLEX(&env);

	return 0;

}

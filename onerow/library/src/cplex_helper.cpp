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

#include <cassert>
#include <cstdio>
#include <iostream>
#include <map>
#include <stdexcept>
#include <algorithm>
#include <ilcplex/cplex.h>
#include <qxx/rational.hpp>
#include <string>
#include <time.h>

#include <onerow/single_row_cut_generator.hpp>
#include <onerow/cplex_helper.hpp>
#include <onerow/geometry.hpp>
#include <onerow/stats.hpp>
#include <onerow/params.hpp>
#include <cstring>

using std::cout;
using std::endl;
using std::string;

static bool debug = false;

CplexHelper::CplexHelper(CPXENVptr _env, CPXLPptr _lp) :
		env(_env), lp(_lp), is_integer(0), n_cuts(0), n_rows(0), ub(0), lb(0),
		cstat(0), cplex_rows(0), first_solution(0), current_solution(0),
		optimal_solution(0), current_round(0), n_good_rows(-1), reduced_costs(0)
{
}


CplexHelper::~CplexHelper()
{
	if (cplex_rows)
	{
		for (int i = 0; i < n_rows; i++)
		{
			delete[] cplex_rows[i].pi;
			delete[] cplex_rows[i].indices;
		}

		delete[] cplex_rows;
	}

	if(ub) delete[] ub;
	if(lb) delete[] lb;
	if(cstat) delete[] cstat;

	if (is_integer)
		delete[] is_integer;
}


#define return_if_neq(a,b) if((a)<(b)) return true; if((a)>(b)) return false;

bool CplexRow::operator<(const CplexRow &other) const
{
	return_if_neq(fabs(pi_zero-0.5), fabs(other.pi_zero-0.5));

	return_if_neq(nz, other.nz);

	for (int i = 0; i < nz; i++)
	{
		return_if_neq(indices[i], other.indices[i]);
		double diff = pi[i] - other.pi[i];
		if(diff < -ZERO_CUTOFF) return true;
	}

	return false;
}

void CplexRow::print(double *solution)
{
	for (int i = 0; i < nz; i++)
	{
		printf("%lf x%d ", pi[i], indices[i]);

		if (solution)
			printf("(%lf) ", solution[indices[i]]);
	}

	printf("<= %lf ", pi_zero);

	if (solution)
		printf("violation=%lf", get_violation(solution));

	printf("\n");
}

double CplexRow::get_violation(double *solution)
{
	double v = 0;
	for (int i = 0; i < nz; i++)
	{
		v += pi[i] * solution[indices[i]];
	}
	v -= pi_zero;

	return v;
}


Constraint CplexHelper::cplex_row_to_constraint(const CplexRow &cplex_row)
{
	Constraint constraint;
	constraint.pi.resize(n_cols);

	rational pi_zero = rational(cplex_row.pi_zero);

	for (int j = 0; j < cplex_row.nz; j++)
	{
		int index = cplex_row.indices[j];
		rational pij = rational(cplex_row.pi[j]);

		if (cstat[index] == CPX_AT_LOWER)
			pi_zero -= rational(lb[index]) * pij;

		if (cstat[index] == CPX_AT_UPPER)
		{
			pi_zero -= rational(ub[index]) * pij;
			pij = -pij;
		}

		constraint.pi.push_nz(cplex_row.indices[j], pij.reduce(REDUCE_FACTOR_COEFFICIENT));
	}

	constraint.pi_zero = rational(pi_zero).reduce(REDUCE_FACTOR_RHS);

	return constraint;
}


CplexRow CplexHelper::constraint_to_cplex_row(const Constraint &cut)
{
	CplexRow cplex_row;
	cplex_row.pi = new double[n_cols];
	cplex_row.indices = new int[cut.pi.nz()];
	cplex_row.pi_zero = cut.pi_zero.get_double();
	cplex_row.depth = cut.depth;
	cplex_row.nz = cut.pi.nz();

	double max_pi = -INFINITY;
	double min_pi = INFINITY;

	int cut_nz = cut.pi.nz();
	for (int j = 0; j < cut_nz; j++)
	{
		int index = cut.pi.index(j);
		double pij = cut.pi.value(j).get_double();

		if (fabs(pij) < ZERO_CUTOFF)
			pij = 0;

		max_pi = std::max(max_pi, fabs(pij));

		if (fabs(pij) > 0)
			min_pi = std::min(min_pi, fabs(pij));

		if (cstat[index] == CPX_AT_LOWER)
			cplex_row.pi_zero += lb[index] * pij;

		if (cstat[index] == CPX_AT_UPPER)
		{
			pij = -pij;
			cplex_row.pi_zero += ub[index] * pij;
		}

		cplex_row.indices[j] = index;
		cplex_row.pi[j] = pij;
	}

	cplex_row.dynamism = max_pi / min_pi;

	return cplex_row;
}


void CplexHelper::add_cut(Constraint *cut)
{

	CplexRow cplex_row = constraint_to_cplex_row(*cut);
	double violation = cplex_row.get_violation(current_solution);

	if (optimal_solution)
		assert(cplex_row.get_violation(optimal_solution) <= MIN_CUT_VIOLATION);

	if (first_solution)
		assert(cplex_row.get_violation(first_solution) >= MIN_CUT_VIOLATION);

	#pragma omp critical
	if (cplex_row.dynamism < MAX_CUT_DYNAMISM && violation >= MIN_CUT_VIOLATION)
	{
		auto p = cut_buffer.insert(cplex_row);

		// duplicate cut
		if (!p.second)
		{
			delete[] cplex_row.pi;
			delete[] cplex_row.indices;
		}

		if (cut_buffer.size() >= MAX_CUT_BUFFER_SIZE)
		{
			flush_cuts();
			solve(false);
		}
	}
	// rejected cut
	else
	{
		delete[] cplex_row.pi;
		delete[] cplex_row.indices;
	}

	delete cut;
}


void CplexHelper::flush_cuts()
{
	int begin = 0;
	char sense = 'L';

	for (CplexRow cplex_row : cut_buffer)
	{
		Stats::add_cut(cplex_row.depth);

		total_cuts++;

		CPXaddrows(env, lp, 0, 1, cplex_row.nz, &cplex_row.pi_zero, &sense,
				&begin, cplex_row.indices, cplex_row.pi, NULL, NULL);

		delete[] cplex_row.pi;
		delete[] cplex_row.indices;
	}

	cut_buffer.clear();
}


Row* CplexHelper::get_tableau_row(int index)
{
	Row *row = new Row;
	row->basic_var_index = cplex_rows[index].head;
	row->is_integer = is_integer;
	row->c = cplex_row_to_constraint(cplex_rows[index]);

	row->reduced_costs = reduced_costs;
	row->cost_cutoff = cost_cutoff;

	if (optimal_solution)
		assert(cplex_rows[index].get_violation(optimal_solution) <= 0.001);

	if (debug)
	{
		printf("ROW %d\n", index);
		cplex_rows[index].print(optimal_solution);
	}

	return row;
}

void CplexHelper::solve(bool should_end_round)
{
	// Optimize
	int status = CPXlpopt(env, lp);
	if (status)
	{
		fprintf(stderr, "Could not optimize (%d)\n", status);
		throw std::runtime_error("CplexHelper::solve_mip");
	}

	// Get status
	char buffer[512];
	double objval;

	status = CPXgetstat(env, lp);

	CPXgetstatstring(env, status, buffer);
	CPXgetobjval(env, lp, &objval);
	time_printf("	 %.6lf [%s]					   \n", objval, buffer);

	assert(status == 1);

	// Store current solution
	if (!current_solution)
		current_solution = new double[n_cols];

	CPXgetx(env, lp, current_solution, 0, n_cols - 1);

	// During first round, store basis and fractional solution
	if (current_round == 0)
	{
		read_basis();

		first_solution = new double[n_cols];
		for (int i = 0; i < n_cols; i++)
			first_solution[i] = current_solution[i];
	}

	if (should_end_round)
		Stats::set_solution(current_round++, objval, string(buffer));
}


void CplexHelper::read_columns()
{
	n_rows = CPXgetnumrows(env, lp);
	n_cols = CPXgetnumcols(env, lp);

	is_integer = new bool[n_cols];

	char ctype[n_cols];
	CPXgetctype(env, lp, ctype, 0, n_cols - 1);

	for (int i = 0; i < n_cols; i++)
	{
		switch (ctype[i])
		{
		case CPX_BINARY:
		case CPX_INTEGER:
		case CPX_SEMIINT:
			is_integer[i] = true;
			break;

		default:
			is_integer[i] = false;
			break;
		}
	}

	time_printf("Fetched %d rows, %d cols.\n", n_rows, n_cols);
}


void CplexHelper::read_basis()
{
	time_printf("Reading basis...\n");

	ub = new double[n_cols];
	lb = new double[n_cols];
	cstat = new int[n_cols];

	int *head = new int[n_rows];
	int *rstat = new int[n_rows];
	double *rhs = new double[n_rows];

	CPXgetbhead(env, lp, head, rhs);
	CPXgetub(env, lp, ub, 0, n_cols - 1);
	CPXgetlb(env, lp, lb, 0, n_cols - 1);
	CPXgetbase(env, lp, cstat, rstat);

	reduced_costs = new double[n_cols];
	CPXgetdj(env, lp, reduced_costs, 0, n_cols-1);

	cost_cutoff = -INFINITY;
	double *costs_copy = new double[n_cols];
	memcpy(costs_copy, reduced_costs, sizeof(double) * n_cols);
	std::sort(costs_copy, costs_copy + n_cols, std::greater<double>());
	if(n_cols > MAX_R1_RAYS) cost_cutoff = costs_copy[MAX_R1_RAYS];
	delete costs_copy;

	cplex_rows = new CplexRow[n_rows];
	assert(cplex_rows != 0);

	eta_reset();
	eta_count = 0;
	eta_total = n_rows;
	std::thread eta(&CplexHelper::eta_print, this);

	for (int i = 0; i < n_rows; i++)
	{
		int nz = 0;
		double pi[n_cols];

		CPXbinvarow(env, lp, i, pi);

		for (int j = 0; j < n_cols; j++)
		{
			if (fabs(pi[j]) < ZERO_CUTOFF)
				continue;

			if (cstat[j] == CPX_AT_LOWER)
				rhs[i] += lb[j] * pi[j];

			if (cstat[j] == CPX_AT_UPPER)
				rhs[i] += ub[j] * pi[j];

			nz++;
		}

		cplex_rows[i].nz = nz;
		cplex_rows[i].depth = 0;
		cplex_rows[i].pi = new double[nz];
		cplex_rows[i].indices = new int[nz];
		cplex_rows[i].pi_zero = rhs[i];
		cplex_rows[i].head = head[i];

		if(fabs(cplex_rows[i].pi_zero) < ZERO_CUTOFF)
			cplex_rows[i].pi_zero = 0;

		int k = 0;
		for (int j = 0; j < n_cols; j++)
		{
			if (fabs(pi[j]) < ZERO_CUTOFF)
				continue;
			cplex_rows[i].pi[k] = pi[j];
			cplex_rows[i].indices[k++] = j;
		}

		eta_count++;
	}

	eta.join();

	delete[] head;
	delete[] rstat;
	delete[] rhs;
}


void CplexHelper::print_basis()
{
	int n_rows = CPXgetnumrows(env, lp);
	int n_cols = CPXgetnumcols(env, lp);

	double y[n_rows];
	double ub[n_cols];
	double lb[n_cols];
	double dj[n_cols];

	CPXgetub(env, lp, ub, 0, n_cols - 1);
	CPXgetlb(env, lp, lb, 0, n_cols - 1);
	CPXgetdj(env, lp, dj, 0, n_cols - 1);

	cout << "basis inverse:" << endl;
	for (int i = 0; i < n_rows; i++)
	{
		CPXbinvrow(env, lp, i, y);
		for (int k = 0; k < n_rows; k++)
			cout << rational(y[k]) << " ";
		cout << endl;
	}

	int cstat[n_cols];
	int rstat[n_rows];
	CPXgetbase(env, lp, cstat, rstat);

	cout << "column status:" << endl;
	for (int i = 0; i < n_cols; i++)
	{
		cout << i << ": " << cstat[i] << " " << lb[i] << "..." << ub[i] << " "
				<< dj[i] << " " << (is_integer[i] ? "int" : "cont") << endl;
	}

}


void CplexHelper::print_solution(double *x)
{
	for (int i = 0; i < n_cols; i++)
		if (fabs(x[i]) > ZERO_CUTOFF)
			time_printf("	 x%d = %.6lf\n", i, x[i]);
}

void CplexHelper::eta_reset()
{
	time(&eta_start);
}

void CplexHelper::eta_print()
{
	while (true)
	{
		for(int i = 0; i < ETA_UPDATE_INTERVAL; i++)
		{
			if(eta_count >= eta_total) goto FINISHED;
			sleep(1);
		}

		if (eta_count == 0)
		{
			time_printf("%3.0f%%  ETA: unknown\n", 0.0);
			continue;
		}

		if (eta_count >= eta_total)
			break;

		time_t eta_now;
		time(&eta_now);

		double diff = difftime(eta_now, eta_start);
		double eta = diff / eta_count * eta_total;

		time_t eta_date = eta_start + eta;
		tm *ttm = localtime(&eta_date);

		time_printf("%3.0f%%  ETA: %04d-%02d-%02d %02d:%02d:%02d  %d / %d\n",
					100.0 * eta_count / eta_total,
					ttm->tm_year + 1900, ttm->tm_mon + 1, ttm->tm_mday,
					ttm->tm_hour, ttm->tm_min, ttm->tm_sec,
					eta_count, eta_total);
	}
FINISHED:;
}

void CplexHelper::find_good_rows(int max_rows)
{
	bool *is_good = new bool[n_rows];
	double *fractionality = new double[n_rows];

	time_printf("Finding interesting rows...\n");

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < n_rows; i++)
	{
		Row *row = get_tableau_row(i);

		fractionality[i] = row->c.pi_zero.frac().get_double();
		fractionality[i] = fabs(fractionality[i] - 0.5);

		is_good[i] = true;
		if (row->c.pi_zero.frac() == 0) is_good[i] = false;
		if (!row->is_integer[row->basic_var_index]) is_good[i] = false;

		delete row;
	}

	if(max_rows > 0)
	{
		double frac_cutoff = 1.0;
		std::sort(fractionality, fractionality + n_rows);
		if (n_rows > max_rows) frac_cutoff = fractionality[max_rows];

		for (int i = 0; i < n_rows; i++)
			if (fractionality[i] > frac_cutoff)
				is_good[i] = false;
	}

	n_good_rows = 0;
	good_rows = new int[n_rows];
	for (int i = 0; i < n_rows; i++)
	{
		if (!is_good[i]) continue;
		good_rows[n_good_rows++] = i;
		if(max_rows > 0 && n_good_rows >= max_rows) break;
	}

	delete is_good;
	delete fractionality;
	time_printf("	 %d rows found\n", n_good_rows);
}

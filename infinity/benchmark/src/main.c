/* Copyright (c) 2015 Alinson Xavier
 *
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

#define _XOPEN_SOURCE 500

#include <getopt.h>
#include <unistd.h>
#include <time.h>

#include <multirow/cg.h>
#include <multirow/mir.h>
#include <multirow/stats.h>
#include <multirow/util.h>

#include <infinity/greedy.h>

int ENABLE_LIFTING = 0;

int MIN_N_ROWS = 2;
int MAX_N_ROWS = 2;

int DUMP_CUT = 0;
int DUMP_CUT_N = 0;

int GENERATE_MIR = 0;
int GENERATE_GREEDY = 0;
int KEEP_INTEGRALITY = 0;
int N_ROUNDS = 1;

char BASIS_FILENAME[1000] = {0};
char PROBLEM_FILENAME[1000] = {0};
char KNOWN_SOLUTION_FILENAME[1000] = {0};

char OUTPUT_SOLUTION_FILENAME[1000] = {0};
char OUTPUT_BASIS_FILENAME[1000] = {0};
char LOG_FILENAME[1000] = {0};
char STATS_FILENAME[1000] = {0};

FILE *LOG_FILE;

int BOOST_VAR = -1;
double BOOST_FACTOR = 0.01;

#define OPTION_WRITE_BASIS 1000
#define OPTION_WRITE_SOLUTION 1001

static const struct option options_tab[] =
{
    {"help", no_argument, 0, 'h'},
    {"problem", required_argument, 0, 'p'},
    {"solution", required_argument, 0, 'x'},
    {"mir", no_argument, 0, 'm'},
    {"greedy", no_argument, 0, 'g'},
    {"rounds", required_argument, 0, 'r'},
    {"keep-integrality", no_argument, 0, 'k'},
    {"write-solution", required_argument, 0, OPTION_WRITE_SOLUTION},
    {"write-basis", required_argument, 0, OPTION_WRITE_BASIS},
    {"basis", required_argument, 0, 'b'},
    {"log", required_argument, 0, 'l'},
    {"stats", required_argument, 0, 's'},
    {"boost", required_argument, 0, 't'},
    {"lift", no_argument, 0, 'i'},
    {"rows", required_argument, 0, 'a'},
    {0, 0, 0, 0}
};

static void print_usage(char **argv)
{
    printf("Usage: %s [OPTION]...\n", argv[0]);
    printf("Solves the given MILP using specified cutting-planes.\n\n");

    printf("Parameters:\n");
    printf("%4s %-20s %s\n", "-b", "--basis=FILE",
           "BAS file containing an optimal basis for the linear relaxation of "
           "the problem");
    printf("%4s %-20s %s\n", "-g", "--greedy",
           "generate greedy intersection cuts");
    printf("%4s %-20s %s\n", "-k", "--keep-integrality",
           "do not relax integrality of variables");
    printf("%4s %-20s %s\n", "-l", "--log=FILE",
           "write log to the specified file");
    printf("%4s %-20s %s\n", "-m", "--mir", "generate MIR cuts");
    printf("%4s %-20s %s\n", "-p", "--problem=FILE", "problem to be solved");
    printf("%4s %-20s %s\n", "-r", "--rounds=NUM",
           "number of rounds of cutting planes to add");
    printf("%4s %-20s %s\n", "-s", "--stats=FILE",
           "write statistics to the specified file");
    printf("%4s %-20s %s\n", "-x", "--solution=FILE",
           "known integral solution (to check the validity of the cuts)");
    printf("%4s %-20s %s\n", "", "--write-solution=FILE",
           "write solution found at the end of the procedure to given file");
    printf("%4s %-20s %s\n", "", "--write-basis=FILE",
           "write optimal LP basis to given file");
    printf("%4s %-20s %s\n", "", "--lift", "enable trivial lifting");
    printf("%4s %-20s %s\n", "", "--rows=N",
            "generate multi-row cuts from up to N rows");
}

static int parse_args(int argc,
                      char **argv)
{
    int rval = 0;
    opterr = 0;

    while (1)
    {
        int c = 0;
        int option_index = 0;
        c = getopt_long(argc, argv, "l:p:s:x:hmgr:kb:t:ia:", options_tab,
                        &option_index);

        if (c < 0) break;

        switch (c)
        {
        case 'l':
            strcpy(LOG_FILENAME, optarg);
            break;

        case 'b':
            strcpy(BASIS_FILENAME, optarg);
            break;

        case 'g':
            GENERATE_GREEDY = 1;
            break;

        case 'k':
            KEEP_INTEGRALITY = 1;
            break;

        case 'm':
            GENERATE_MIR = 1;
            break;

        case 't':
            BOOST_VAR = atoi(optarg);
            break;

        case 'r':
            N_ROUNDS = atoi(optarg);
            break;

        case 'a':
            MAX_N_ROWS = atoi(optarg);
            break;

        case OPTION_WRITE_SOLUTION:
            strcpy(OUTPUT_SOLUTION_FILENAME, optarg);
            break;

        case OPTION_WRITE_BASIS:
            strcpy(OUTPUT_BASIS_FILENAME, optarg);
            break;

        case 'p':
            strcpy(PROBLEM_FILENAME, optarg);
            break;

        case 's':
            strcpy(STATS_FILENAME, optarg);
            break;

        case 'x':
            strcpy(KNOWN_SOLUTION_FILENAME, optarg);
            break;

        case 'h':
            print_usage(argv);
            exit(0);

        case 'i':
            ENABLE_LIFTING = 1;
            break;

        case ':':
            fprintf(stderr, "%s: option '-%c' requires an argument\n",
                    argv[0], optopt);
            rval = 1;
            goto CLEANUP;

        case '?':
        default:
            fprintf(stderr, "%s: option '-%c' is invalid\n", argv[0],
                    optopt);
            rval = 1;
            goto CLEANUP;

        }
    }
    if ((strlen(PROBLEM_FILENAME) == 0))
    {
        fprintf(stderr, "You must specify the problem.\n");
        rval = 1;
    }

    if (KEEP_INTEGRALITY && (GENERATE_GREEDY + GENERATE_MIR > 0))
    {
        fprintf(stderr, "Cutting planes cannot be added when integrality is "
                "kept\n");
        rval = 1;
    }

    if (N_ROUNDS < 1)
    {
        fprintf(stderr, "Invalid number of rounds.\n");
        rval = 1;
    }

CLEANUP:
    if (rval)
        fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
    return rval;
}

void print_header(int argc,
                  char *const *argv)
{
    char hostname[3000];
    gethostname(hostname, 1024);

    time_t now;
    time(&now);
    struct tm *ttm = localtime(&now);

    time_printf("multirow\n");
    time_printf("%s %04d-%02d-%02d %02d:%02d\n", hostname, ttm->tm_year + 1900,
                ttm->tm_mon + 1, ttm->tm_mday, ttm->tm_hour, ttm->tm_min);

    time_printf("Compile-time parameters:\n");
    time_printf("    EPSILON: %e\n", EPSILON);
    time_printf("    GREEDY_BIG_E: %e\n", GREEDY_BIG_E);
    time_printf("    GREEDY_MAX_GAP: %e\n", GREEDY_MAX_GAP);

    char cmdline[5000] = {0};
    for (int i = 0; i < argc; i++)
        sprintf(cmdline + strlen(cmdline), "%s ", argv[i]);

    time_printf("Command line arguments:\n");
    time_printf("    %s\n", cmdline);
}

int main(int argc,
         char **argv)
{
    int rval = 0;
    double *x = 0;
    struct CG *cg = 0;

    struct LP lp;
    char *column_types = 0;

    rval = parse_args(argc, argv);
    if (rval) return 1;

    if (LOG_FILENAME[0])
    {
        LOG_FILE = fopen(LOG_FILENAME, "w");
        abort_if(!LOG_FILE, "could not open log file");
    }

    if_info_level
    {
        print_header(argc, argv);
    }

    STATS_init();
    STATS_set_input_filename(PROBLEM_FILENAME);
    progress_title(PROBLEM_FILENAME);

    rval = LP_open(&lp);
    abort_if(rval, "LP_open failed");

    rval = LP_create(&lp, "multirow");
    abort_if(rval, "LP_create failed");

    rval = LP_read_problem(&lp, PROBLEM_FILENAME);
    abort_if(rval, "LP_read_problem failed");

    int ncols = LP_get_num_cols(&lp);
    int nrows = LP_get_num_rows(&lp);

    log_info("    %d rows, %d cols\n", nrows, ncols);

    x = (double *) malloc(ncols * sizeof(double));
    abort_if(!x, "could not allocate x");

    if (!KEEP_INTEGRALITY)
    {
        column_types = (char *) malloc(ncols * sizeof(char));
        abort_if(!column_types, "could not allocate column_types");

        log_info("Storing column types...\n");
        rval = LP_get_column_types(&lp, column_types);
        abort_if(rval, "LP_get_column_types failed");

        log_info("Relaxing integrality...\n");
        rval = LP_relax(&lp);
        abort_if(rval, "LP_relax failed");

        log_info("Disabling presolve...\n");
        LP_disable_presolve(&lp);
    }

    if(BASIS_FILENAME[0])
    {
        rval = LP_read_basis(&lp, BASIS_FILENAME);
        abort_if(rval, "LP_read_basis failed");
    }

    log_info("Optimizing...\n");
    int infeasible;
    rval = LP_optimize(&lp, &infeasible);
    abort_if(rval, "LP_optimize failed");

    double cost = 0, xboost = 0;
    if(BOOST_VAR > 0)
    {
        rval = LP_get_x(&lp, x);
        abort_if(rval, "LP_get_x failed");
        xboost = x[BOOST_VAR];

        for (int i = 0; i < ncols; i++)
            log_info("x[%3d] = %12.8lf\n", i, x[i]);
    }


    double obj;
    rval = LP_get_obj_val(&lp, &obj);
    abort_if(rval, "LP_get_obj_val failed");

    log_info("    opt = %lf\n", obj);

    STATS_set_obj_value(obj);
    STATS_finish_round();

    if(BOOST_VAR >= 0)
        log_info("Boosting variable %d by %.2lf\n", BOOST_VAR, BOOST_FACTOR);

    if(OUTPUT_BASIS_FILENAME[0])
    {
        rval = LP_write_basis(&lp, OUTPUT_BASIS_FILENAME);
        abort_if(rval, "LP_write_basis failed");
    }

    if(GENERATE_MIR || GENERATE_GREEDY)
    {
        cg = (struct CG *) malloc(sizeof(struct CG));
        abort_if(!cg, "could not allocate cg");

        log_info("Reading tableau rows...\n");
        rval = CG_init(&lp, column_types, cg);
        abort_if(rval, "CG_init failed");

        if (strlen(KNOWN_SOLUTION_FILENAME) > 0)
        {
            if(access(KNOWN_SOLUTION_FILENAME, F_OK) != -1)
            {
                rval = LP_read_solution(&lp, KNOWN_SOLUTION_FILENAME, x);
                abort_if(rval, "LP_read_solution failed");

                rval = CG_set_integral_solution(cg, x);
                abort_if(rval, "CG_set_integral_solution failed");
            }
            else
            {
                log_error("ERROR: Could not read solution file!\n");
            }
        }

        rval = LP_get_x(&lp, x);
        abort_if(rval, "LP_get_x failed");

        rval = CG_set_basic_solution(cg, x);
        abort_if(rval, "CG_set_basic_solution failed");

        if (GENERATE_MIR)
        {
            log_info("Adding MIR cuts...\n");
            rval = CG_add_single_row_cuts(
                cg,
                (SingleRowGeneratorCallback)
                MIR_generate_cut);
            abort_if(rval, "CG_add_single_row_cuts failed");

            log_info("Optimizing...\n");
            rval = LP_optimize(&lp, &infeasible);
            abort_if(rval, "LP_optimize failed");

            rval = LP_get_obj_val(&lp, &obj);
            abort_if(rval, "LP_get_obj_val failed");

            log_info("    opt = %lf\n", obj);
            STATS_set_obj_value(obj);

            STATS_finish_round();
        }

        if (GENERATE_GREEDY)
        {
            for(int k = MIN_N_ROWS; k <= MAX_N_ROWS; k++)
            {
                log_info("Adding greedy intersection cuts (%d rows)...\n", k);

                rval = CG_add_multirow_cuts(cg, k, (MultirowGeneratorCallback)
                                            GREEDY_generate_cut);
                abort_if(rval, "CG_add_multirow_cuts failed");
            }

            log_info("Optimizing...\n");
            rval = LP_optimize(&lp, &infeasible);
            abort_if(rval, "LP_optimize failed");

            rval = LP_get_obj_val(&lp, &obj);
            abort_if(rval, "LP_get_obj_val failed");

            log_info("    opt = %lf\n", obj);
            STATS_set_obj_value(obj);

            STATS_finish_round();
        }

        CG_free(cg);
    }


    if(BOOST_VAR > 0)
    {
        FILE *boost = fopen("boost.csv", "a");
        fprintf(boost, "%d, %.8lf, %.8lf, %.8lf\n", BOOST_VAR, obj, cost,
                xboost);
    }

    if(OUTPUT_SOLUTION_FILENAME[0])
    {
        rval = LP_write_solution(&lp, OUTPUT_SOLUTION_FILENAME);
        abort_if(rval, "LP_write_solution failed");
    }

    if(STATS_FILENAME[0])
    {
        log_info("Writing stats to file %s...\n", STATS_FILENAME);
        rval = STATS_print_yaml(STATS_FILENAME);
        abort_if(rval, "STATS_print_yaml failed");
    }

CLEANUP:
    if (LOG_FILE) fclose(LOG_FILE);
    if (x) free(x);
    if (column_types) free(column_types);
    LP_free(&lp);

    return rval;
}

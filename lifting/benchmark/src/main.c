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

#include <stdio.h>
#include <stdarg.h>
#include <getopt.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>

#include <multirow/util.h>
#include <multirow/double.h>
#include <multirow/lfree2d.h>
#include <lifting/lifting.h>
#include <lifting/lifting-mip.h>

char LOG_FILENAME[1000] = {0};
char STATS_FILENAME[1000] = {0};
char SETS_FILENAME[1000] = {0};
char ANSWERS_FILENAME[1000] = {0};
unsigned int SEED = 0;

#define ALGORITHM_BOUND 0
#define ALGORITHM_NAIVE 1
#define ALGORITHM_MIP 2

int SELECT_NAIVE_ALGORITHM = 0;
int SELECT_BOUND_ALGORITHM = 0;
int SELECT_MIP_ALGORITHM = 0;

int ENABLE_PREPROCESSING = 0;
int ENABLE_SHEAR = 0;

int CHECK_ANSWERS = 0;
int WRITE_ANSWERS = 0;

int USE_FIXED_BOUNDS = 0;
int NAIVE_BIG_M = 0;
int N_SAMPLES_PER_SET = 10;

int N_ANSWERS = 0;
double ANSWERS[100000];

double PRE_M[100000];
double CENTER[100000];

FILE *LOG_FILE;
FILE *STATS_FILE;
FILE *ANSWERS_FILE;

int BOUNDING_BOX_PADDING = 5;

static const struct option options_tab[] =
{
    {"help", no_argument, 0, 'h'},
    {"sets", required_argument, 0, 'b'},
    {"log", required_argument, 0, 'l'},
    {"stats", required_argument, 0, 'o'},
    {"seed", required_argument, 0, 's'},
    {"fixed-bounds", required_argument, 0, 'f'},
    {"naive", no_argument, 0, 'n'},
    {"bound", no_argument, 0, 'u'},
    {"preprocess", no_argument, 0, 'p'},
    {"shear", no_argument, 0, 'e'},
    {"write-answers", required_argument, 0, 'w'},
    {"check-answers", required_argument, 0, 'c'},
    {"samples", required_argument, 0, 'a'},
    {"mip", no_argument, 0, 'm'},
    {0, 0, 0, 0}
};

static void print_usage(char **argv)
{
    printf("Usage: %s [OPTION]...\n", argv[0]);
    printf("Performs trivial lifting using different algorithms.\n\n");

    printf("Parameters:\n");
    printf("%4s %-20s %s\n", "-l", "--log=FILE",
           "write log to the specified file");
    printf("%4s %-20s %s\n", "-b", "--sets=FILE",
            "file containing lattice-free sets");
    printf("%4s %-20s %s\n", "-o", "--stats=FILE",
           "write statistics to the specified file");
    printf("%4s %-20s %s\n", "-s", "--seed=NUM",
            "random seed");
    printf("%4s %-20s %s\n", "-f", "--fixed-bounds=NUM",
            "user fixed bounds for naive algorithm");
    printf("%4s %-20s %s\n", "-n", "--naive",
            "select naive algorithm");
    printf("%4s %-20s %s\n", "-m", "--mip",
            "select MIP algorithm");
    printf("%4s %-20s %s\n", "-u", "--bound",
            "select bound algorithm");
    printf("%4s %-20s %s\n", "-p", "--preprocess",
            "enable pre-processing step in bound algorithm");
    printf("%4s %-20s %s\n", "-e", "--shear",
            "apply shear transformation to sets");
    printf("%4s %-20s %s\n", "-w", "--write-answers=FILE",
            "write computed coefficients to given file");
    printf("%4s %-20s %s\n", "-c", "--check-answers=FILE",
            "check computed coefficients against given file");
    printf("%4s %-20s %s\n", "-a", "--samples=NUM",
            "use specified number of samples per set");
}

void stats_printf(const char *fmt,
                 ...)
{
    if(STATS_FILE)
    {
        va_list args;
        va_start(args, fmt);
        vfprintf(STATS_FILE, fmt, args);
        va_end(args);

        fflush(STATS_FILE);
    }
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
        c = getopt_long(argc, argv, "hb:k:s:f:o:nupew:c:a:m", options_tab,
                        &option_index);

        if (c < 0) break;

        switch (c)
        {
        case 'l':
            strcpy(LOG_FILENAME, optarg);
            break;

        case 's':
        {
            int count = sscanf(optarg, "%u", &SEED);
            abort_if(count != 1, "invalid seed");
            break;
        }

        case 'a':
        {
            int count = sscanf(optarg, "%d", &N_SAMPLES_PER_SET);
            abort_if(count != 1, "invalid number of samples");
            abort_if(N_SAMPLES_PER_SET <= 0, "invalid number of samples");
            break;
        }

        case 'f':
        {
            USE_FIXED_BOUNDS = 1;
            int count = sscanf(optarg, "%d", &NAIVE_BIG_M);
            abort_if(count != 1, "invalid fixed bound");
            break;
        }

        case 'o':
            strcpy(STATS_FILENAME, optarg);
            break;

        case 'b':
            strcpy(SETS_FILENAME, optarg);
            break;

        case 'w':
            strcpy(ANSWERS_FILENAME, optarg);
            WRITE_ANSWERS = 1;
            break;

        case 'c':
            strcpy(ANSWERS_FILENAME, optarg);
            CHECK_ANSWERS = 1;
            break;

        case 'n':
            SELECT_NAIVE_ALGORITHM = 1;
            break;

        case 'm':
            SELECT_MIP_ALGORITHM = 1;
            break;

        case 'u':
            SELECT_BOUND_ALGORITHM = 1;
            break;

        case 'p':
            ENABLE_PREPROCESSING = 1;
            break;

        case 'e':
            ENABLE_SHEAR = 1;
            break;

        case 'h':
            print_usage(argv);
            exit(0);

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

    if ((strlen(SETS_FILENAME) == 0))
    {
        fprintf(stderr, "You must specify a file containing the lattice-free "
                "sets.\n");
        rval = 1;
    }

    if (SELECT_NAIVE_ALGORITHM + SELECT_BOUND_ALGORITHM + SELECT_MIP_ALGORITHM != 1)
    {
        fprintf(stderr, "You must select exactly one algorithm.\n");
        rval = 1;
    }

    if (CHECK_ANSWERS + WRITE_ANSWERS > 1)
    {
        fprintf(stderr, "Cannot write and check answers at same time.\n");
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

    time_printf("benchmark.run\n");
    time_printf("%s %04d-%02d-%02d %02d:%02d\n", hostname, ttm->tm_year + 1900,
                ttm->tm_mon + 1, ttm->tm_mday, ttm->tm_hour, ttm->tm_min);

    time_printf("Compile-time parameters:\n");
    time_printf("    EPSILON: %e\n", EPSILON);

    char cmdline[5000] = {0};
    for (int i = 0; i < argc; i++)
        sprintf(cmdline + strlen(cmdline), "%s ", argv[i]);

    time_printf("Command line arguments:\n");
    time_printf("    %s\n", cmdline);
}

int benchmark_set_sample(int algorithm,
                         int set_idx,
                         const struct LFreeSet2D *set,
                         const double *rays,
                         const double *pre_m,
                         const double center,
                         int *lb,
                         int *ub,
                         int current_sample,
                         int *wrong_answer)
{
    int rval = 0;

    for (int i = 0; i < N_RAYS; i++)
    {
        double ray[2] = { rays[2 * i], rays[2 * i + 1] };
        double value;

        if(ENABLE_PREPROCESSING)
        {
            double r0 = ray[0], r1 = ray[1];
            ray[0] = pre_m[0] * r0 + pre_m[2] * r1;
            ray[1] = pre_m[1] * r0 + pre_m[3] * r1;

            ray[0] = ray[0] - floor(set->f[0] + ray[0]);
            ray[1] = ray[1] + floor(center + 0.5 - set->f[1] - ray[1]);
        }

        log_debug("    Ray %d (%.6lf,%.6lf)...\n", i, ray[0], ray[1]);

        switch (algorithm)
        {
            case ALGORITHM_BOUND:
                rval = LIFTING_2D_bound(set->n_halfspaces, set->halfspaces, ray,
                        &value);
                abort_if(rval, "LIFTING_2D_bound failed");
                break;

            case ALGORITHM_NAIVE:
                rval = LIFTING_2D_naive(set->n_halfspaces, set->halfspaces, ray,
                        lb, ub, &value);
                abort_if(rval, "LIFTING_2D_naive failed");
                break;

            case ALGORITHM_MIP:
                rval = LIFTING_2D_mip(set->n_halfspaces, set->halfspaces, ray,
                        &value);
                abort_if(rval, "LIFTING_2D_mip failed");
                break;

            default:
                abort_if(1, "Invalid algorithm");
        }

        if(current_sample == 0)
        {
            if(WRITE_ANSWERS)
            {
                abort_iff(!DOUBLE_geq(value, 0),
                        "value should be non-negative (%.8lf)", value);
                fprintf(ANSWERS_FILE, "%d %d %.20lf\n", set_idx, i, value);
            }

            if(CHECK_ANSWERS)
            {
                int count;
                int expected_set_idx, expected_i;
                double expected_value;

                while(1)
                {
                    count = fscanf(ANSWERS_FILE, "%d %d %lf ",
                            &expected_set_idx, &expected_i, &expected_value);

                    abort_if(count != 3, "error reading answer");
                    if(set_idx == expected_set_idx && i == expected_i) break;
                }

                double delta = fabs(value - expected_value);
                if(delta > 1e-3)
                {
                    log_warn("    wrong answer (set=%d ray=%d answer=%.8lf"
                            " expected=%.8lf delta=%.8lf)\n", set_idx,
                            i, value, expected_value, delta);
                    *wrong_answer = 1;

                    LFREE_2D_print_set(set);
                }
            }

            log_verbose("       %4d %12.8lf\n", j, value);
        }
    }

CLEANUP:
    return rval;
}

int benchmark_set(int algorithm,
                  int set_idx,
                  const struct LFreeSet2D *set,
                  const double *rays,
                  const double *pre_m,
                  const double center,
                  int *wrong_answer)
{
    int rval = 0;
    int lb[2], ub[2];

    if(algorithm == ALGORITHM_NAIVE)
    {
        if (USE_FIXED_BOUNDS)
        {
            ub[0] = ub[1] = NAIVE_BIG_M;
            lb[0] = lb[1] = -NAIVE_BIG_M;
        }
        else
        {
            rval = LFREE_2D_get_bounding_box(set, lb, ub);
            abort_if(rval, "LFREE_2D_get_bounding_box failed");

            ub[0] += BOUNDING_BOX_PADDING;
            ub[1] += BOUNDING_BOX_PADDING;
            lb[0] -= BOUNDING_BOX_PADDING;
            lb[0] -= BOUNDING_BOX_PADDING;
        }

        ub[0] = fmin(ub[0], MAX_BOX_SIZE);
        ub[1] = fmin(ub[1], MAX_BOX_SIZE);
        lb[0] = fmax(lb[0], -MAX_BOX_SIZE);
        lb[1] = fmax(lb[1], -MAX_BOX_SIZE);

        if(ub[0] == MAX_BOX_SIZE || ub[1] == MAX_BOX_SIZE || lb[0] ==
                -MAX_BOX_SIZE || lb[1] == -MAX_BOX_SIZE)
        {
            log_info("    bounding box has been clipped\n");
        }
    }

    for (int k = 0; k < N_SAMPLES_PER_SET; k ++)
    {
        rval = benchmark_set_sample(algorithm, set_idx, set, rays, pre_m,
                center, lb, ub, k, wrong_answer);
        abort_if(rval, "benchmark_set_sample failed");
    }

CLEANUP:
    return rval;
}

int benchmark(int n_sets, struct LFreeSet2D *sets, double *rays,
        int algorithm)
{
    int rval = 0;

    log_info("Running benchmark...\n");

    double total_initial_time = get_user_time();
    stats_printf("cpu_time:\n");

    for (int j = 0; j < n_sets; j++)
    {
        log_debug("Set %d...\n", j);

        double set_initial_time = get_user_time();

        struct LFreeSet2D *set = &sets[j];
        double *pre_m = &PRE_M[j * 4];
        double center = CENTER[j];

        rval = LFREE_2D_print_set(set);
        abort_if(rval, "LFREE_2D_print_set failed");

        int wrong_answer = 0;

        rval = benchmark_set(algorithm, j, set, rays, pre_m, center, &wrong_answer);
        abort_if(rval, "benchmark_set failed");

        double set_duration = get_user_time() - set_initial_time;
        double avg = (set_duration / N_SAMPLES_PER_SET) * 1000;

        if(wrong_answer) avg = 1000000;

        stats_printf("  %d: %.8lf\n", j, avg);
        log_info("  %3d: %12.3lf ms\n", j, avg);
    }

    double total_duration = get_user_time() - total_initial_time;

    log_info("    %.3lf ms per set                     \n",
            total_duration / (n_sets * N_SAMPLES_PER_SET) * 1000);

    if(algorithm == ALGORITHM_MIP)
    {
        log_info("    %.3lf s spent on LP_create\n", MIP_TIME_CREATE);
        log_info("    %.3lf s spent on LP_optimize\n", MIP_TIME_OPTIMIZE);
    }

CLEANUP:
    return rval;
}

int main(int argc, char **argv)
{
    int rval = 0;
    double *rays = 0;
    struct LFreeSet2D sets[MAX_N_SETS];

    rval = parse_args(argc, argv);
    if (rval) return 1;

    print_header(argc, argv);

    if (LOG_FILENAME[0])
    {
        LOG_FILE = fopen(LOG_FILENAME, "w");
        abort_if(!LOG_FILE, "could not open log file");
        log_info("Writing log to file: %s\n", LOG_FILENAME);
    }

    if (STATS_FILENAME[0])
    {
        STATS_FILE = fopen(STATS_FILENAME, "w");
        abort_if(!STATS_FILE, "could not open stats file");
        log_info("Writing stats to file: %s\n", STATS_FILENAME);
    }

    if (WRITE_ANSWERS)
    {
        N_SAMPLES_PER_SET = 1;
        ANSWERS_FILE = fopen(ANSWERS_FILENAME, "w");
        abort_if(!ANSWERS_FILE, "could not open answers file");
        log_info("Writing answers to file: %s\n", ANSWERS_FILENAME);
    }

    if (CHECK_ANSWERS)
    {
        ANSWERS_FILE = fopen(ANSWERS_FILENAME, "r");
        abort_if(!ANSWERS_FILE, "could not open answers file");
        log_info("Reading answers from file: %s\n", ANSWERS_FILENAME);
    }

    if(SEED == 0)
    {
        struct timeval t1;
        gettimeofday(&t1, NULL);
        SEED = (unsigned int) (t1.tv_usec * t1.tv_sec) % 10000;
    }

    log_info("Random seed: %u\n", SEED);
    srand(SEED);

    log_info("Generating %d random rays...\n", N_RAYS);

    rays = (double*) malloc(2 * N_RAYS * sizeof(double));
    abort_if(!rays, "could not allocate rays");

    for (int i = 0; i < N_RAYS; i++)
    {
        double *ray = &rays[2 * i];
        ray[0] = DOUBLE_random(0.0, 1.0);
        ray[1] = DOUBLE_random(0.0, 1.0);
    }

    int algorithm = ALGORITHM_BOUND;
    if(SELECT_MIP_ALGORITHM) algorithm = ALGORITHM_MIP;
    else if(SELECT_NAIVE_ALGORITHM) algorithm = ALGORITHM_NAIVE;

    if(algorithm == ALGORITHM_NAIVE)
    {
        log_info("Enabling naive algorithm\n");

        if(USE_FIXED_BOUNDS)
            log_info("Using fixed big M: %d\n", NAIVE_BIG_M);
        else
            log_info("Enabling bounding boxes\n");
    }
    else
    {
        log_info("Enabling bound algorithm\n");
    }

    log_info("Setting %d samples per set\n", N_SAMPLES_PER_SET);

    if(ENABLE_PREPROCESSING)
        log_info("Enabling pre-processing\n");

    log_info("Reading sets from file...\n");
    FILE *sets_file = fopen(SETS_FILENAME, "r");
    abort_iff(!sets_file, "could not read file %s", SETS_FILENAME);

    int line = 0;
    int n_sets = 0;
    while(!feof(sets_file))
    {
        line++;
        struct LFreeSet2D *set = &sets[n_sets];
        LFREE_2D_init(set, 4, 4, 4);

        rval = LFREE_2D_read_next(sets_file, set);
        abort_iff(rval, "LFREE_2D_read_next failed (line %d)", line);

        if(ENABLE_SHEAR)
        {
            double m[4] = { 51.0, 5.0, 10.0, 1.0 };
            rval = LFREE_2D_transform_set(set, m);
            abort_iff(rval, "LFREE_2D_transform_set failed (line %d)", line);
        }

        double dx = -floor(set->f[0]);
        double dy = -floor(set->f[1]);
        rval = LFREE_2D_translate_set(set, dx, dy);
        abort_iff(rval, "LFREE_2D_translate_set failed (line %d)", line);

        if(ENABLE_PREPROCESSING)
        {
            double *pre_m = &PRE_M[n_sets * 4];
            double *center = &CENTER[n_sets];
            rval = LFREE_2D_preprocess(set, pre_m, center);
            abort_iff(rval, "LFREE_2D_preprocess failed (line %d)", line);
        }

        rval = LFREE_2D_compute_halfspaces(set);
        abort_iff(rval, "LFREE_2D_compute_halfspaces failed (line %d)", line);

        rval = LIFTING_2D_verify(set);
        if(rval)
        {
            log_warn("    skipping invalid set on line %d\n", line);
            continue;
        }

        n_sets++;
        if(n_sets >= MAX_N_SETS) break;
    }

    fclose(sets_file);
    log_info("Successfully read %d sets\n", n_sets);

    rval = benchmark(n_sets, sets, rays, algorithm);
    abort_if(rval, "benchmark failed");

    log_info("Done.\n");

CLEANUP:
    if (LOG_FILE) fclose(LOG_FILE);
    if (STATS_FILE) fclose(STATS_FILE);
    if (ANSWERS_FILE) fclose(ANSWERS_FILE);
    if (rays) free(rays);
    return rval;
}

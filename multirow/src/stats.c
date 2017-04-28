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

#define MAX_LENGTH 1000
#define MAX_ROUNDS 100

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <multirow/util.h>

char input_filename[MAX_LENGTH];

unsigned long long generated_cuts_count[MAX_LENGTH];
unsigned long long added_cuts_count[MAX_LENGTH];

int n_rounds = 0;
unsigned long long total_added_cuts = 0;
unsigned long long total_generated_cuts = 0;

double obj_value[MAX_ROUNDS];

double prev_time;
double runtime[MAX_ROUNDS];

void STATS_init()
{
    for(int i = 0; i < MAX_ROUNDS; i++)
    {
        obj_value[i] = NAN;
        generated_cuts_count[i] = 0;
        added_cuts_count[i] = 0;
        runtime[i] = 0;
    }

    prev_time = get_user_time();
}

void STATS_set_input_filename(char *filename)
{
    strcpy(input_filename, filename);
}

void STATS_set_obj_value(double obj)
{
    obj_value[n_rounds] = obj;
}

void STATS_increment_generated_cuts()
{
    generated_cuts_count[n_rounds]++;
    total_generated_cuts++;
}

void STATS_increment_added_cuts()
{
    added_cuts_count[n_rounds]++;
    total_added_cuts++;
}

void STATS_finish_round()
{
    double now = get_user_time();
    runtime[n_rounds] = now - prev_time;
    prev_time = now;

    n_rounds++;
}

int STATS_print_yaml(char *filename)
{
    int rval = 0;

    FILE *f = fopen(filename, "w");
    abort_iff(!f, "could not open file %s", filename);

    fprintf(f, "input-filename:\n  %s\n", input_filename);

    fprintf(f, "obj-value:\n");
    for(int i = 0; i < n_rounds; i++)
        fprintf(f, "  %d: %.6lf\n", i, obj_value[i]);

    fprintf(f, "generated-cuts:\n");
    for(int i = 0; i < n_rounds; i++)
        fprintf(f, "  %d: %lld\n", i, generated_cuts_count[i]);

    fprintf(f, "added-cuts:\n");
    for(int i = 0; i < n_rounds; i++)
        fprintf(f, "  %d: %lld\n", i, added_cuts_count[i]);

    fprintf(f, "user-cpu-time:\n");
    for(int i = 0; i < n_rounds; i++)
        fprintf(f, "  %d: %.3lf\n", i, runtime[i]);

    fprintf(f, "time-per-cut:\n");
    for(int i = 0; i < n_rounds; i++)
    {
        if(generated_cuts_count[i] > 0)
            fprintf(f, "  %d: %.3lf\n", i, runtime[i] / generated_cuts_count[i]);
    }

    fclose(f);

CLEANUP:
    return rval;
}


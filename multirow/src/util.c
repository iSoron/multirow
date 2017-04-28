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

#include <stdio.h>
#include <sys/resource.h>
#include <sys/ioctl.h>
#include <stdarg.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#include <multirow/util.h>

FILE *LOG_FILE = 0;
double INITIAL_TIME = 0;

double get_user_time()
{
    struct rusage ru;

    getrusage(RUSAGE_SELF, &ru);

    return ((double) ru.ru_utime.tv_sec)
           + ((double) ru.ru_utime.tv_usec) / 1000000.0;
}

double get_real_time()
{
    return (double) time (0);
}

void time_printf(const char *fmt,
                 ...)
{
    if (INITIAL_TIME == 0)
        INITIAL_TIME = get_user_time();

    printf("[%10.2lf] ", get_user_time() - INITIAL_TIME);

    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
    fflush(stdout);

    if(LOG_FILE)
    {
        fprintf(LOG_FILE, "[%10.2lf] ", get_user_time() - INITIAL_TIME);

        va_list args;
        va_start(args, fmt);
        vfprintf(LOG_FILE, fmt, args);
        va_end(args);

        fflush(LOG_FILE);
    }
}

double frac(double x)
{
    return x - floor(x);
}

static long eta_count;
static long eta_total;
static time_t eta_start;
static time_t eta_last;
static char eta_title[1000] = {0};

void progress_reset()
{
    eta_count = 0;
    time(&eta_start);
    eta_last = eta_start;
}

void progress_set_total(long total)
{
    eta_total = total;
}

void progress_increment()
{
    eta_count++;
}

void progress_title(char *title)
{
    strncpy(eta_title, title, 1000);
}

void progress_print()
{
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);

    if (eta_count == 0) return;

    time_t eta_now;
    time(&eta_now);

    if(fabs(difftime(eta_now, eta_last)) < 30.0)
        return;

    eta_last = eta_now;

    double diff = difftime(eta_now, eta_start);
    double eta = diff / eta_count * eta_total;

    int length = 10;

    time_t eta_date = eta_start + eta;
    struct tm *ttm = localtime(&eta_date);

    fprintf(stdout,
            "[%-40s] %*s%3.0f%%  ETA: %04d-%02d-%02d %02d:%02d  %*ld / %*ld\n",
            eta_title, 0, "", 100.0 * eta_count / eta_total, ttm->tm_year +
            1900, ttm->tm_mon + 1, ttm->tm_mday, ttm->tm_hour, ttm->tm_min,
            length, eta_count, length, eta_total);
    fflush(stdout);
}

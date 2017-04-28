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
#ifndef MULTIROW_STATS_H
#define MULTIROW_STATS_H

void STATS_init();

void STATS_set_input_filename(char *filename);

void STATS_set_obj_value(double obj);

void STATS_increment_generated_cuts();

void STATS_increment_added_cuts();

int STATS_print_yaml(char *filename);

void STATS_finish_round();

#endif

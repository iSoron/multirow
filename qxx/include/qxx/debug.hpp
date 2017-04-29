/*
    This file is part of qxx -- matrix algebra in exact arithmetic
    Copyright (C) 2013-2014  Laurent Poirrier

    libp is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with pxx.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef QXX_DEBUG_H
#define QXX_DEBUG_H
#include <cassert>
#include <cstdio>

// **************************************************************************
// 
// **************************************************************************
#define Q_RANGE_CHECK(i, lb, ub)	\
	do { \
		if (((i) < (lb)) || ((i) > (ub))) { \
			fprintf(stderr, "%s:%d: in %s(): range check failed " \
				"(%s not in %s..%s | %d not in %d..%d)\n", \
				__FILE__, __LINE__, __func__, \
				#i, #lb, #ub, (i), (lb), (ub)); \
			assert(!"range check"); \
		} \
	} while(0)

#define Q_MATCH(a, b)	\
	do { \
		if ((a) != (b)) { \
			fprintf(stderr, "%s:%d: in %s(): match failed " \
				"(%s != %s | %d != %d)\n", \
				__FILE__, __LINE__, __func__, \
				#a, #b, (a), (b)); \
			assert(!"match"); \
		} \
	} while(0)


// **************************************************************************
// 
// **************************************************************************

#endif



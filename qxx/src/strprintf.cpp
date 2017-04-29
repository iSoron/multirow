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
#include <cstdio>
#include <cstdarg>
#include <gmp.h>
#include "qxx/strprintf.hpp"


namespace q {

// ***************************************************************************
// 
// ***************************************************************************
void format(std::string &s, const char *fmt, ...)
{
	va_list ap;
	int l;
	
	va_start(ap, fmt);
	l = gmp_vsnprintf(NULL, 0, fmt, ap);
	va_end(ap);
	
	if (l <= 0) {
		s = "";
		return;
	}
	
	char *store = new char[l + 1];
	
	va_start(ap, fmt);
	gmp_vsnprintf(store, l + 1, fmt, ap);
	va_end(ap);
	
	s = store;
	
	delete[] store;
}


// ***************************************************************************
// 
// ***************************************************************************
std::string strprintf(const char *fmt, ...)
{
	va_list ap;
	int l;
	
	va_start(ap, fmt);
	l = gmp_vsnprintf(NULL, 0, fmt, ap);
	va_end(ap);
	
	if (l <= 0)
		return(std::string());
	
	char *store = new char[l + 1];
	
	va_start(ap, fmt);
	gmp_vsnprintf(store, l + 1, fmt, ap);
	va_end(ap);
	
	std::string s = store;
	
	delete[] store;
	
	return(s);
}



} // namespace


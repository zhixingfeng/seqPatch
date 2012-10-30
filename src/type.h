#ifndef TYPE_H
#define TYPE_H

#include "stl.h"

#if defined(_MSC_VER) || defined(__BORLANDC__)
	#define WINDOWS
#endif

#ifdef WINDOWS
	typedef __int64 int64;
	typedef unsigned __int64 uint64;
#else
	typedef long long int64;
	typedef unsigned long long uint64;
#endif

typedef unsigned char uchar;
typedef unsigned int uint;
		
#define __DEBUG

#ifdef __ASSERT
	#undef __ASSERT
#endif

#ifdef __DEBUG
	#define __ASSERT(b, msg) if (!(b)) panic(msg);
#else
	#define __ASSERT
#endif

inline void panic(const string errormsg){
	cout << errormsg.c_str() << endl;
	exit(1);
}

#endif //TYPE_H

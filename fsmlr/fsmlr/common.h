#ifndef _COMMON_H_
#define _COMMON_H_

#define INT8 signed char
#define UINT8 unsigned char
#define INT16 signed short int
#define UINT16 unsigned short int
#define INT32 signed long int
#define UINT32 unsigned long int
#ifdef WIN32
	#define INT64 __int64
#else
	#define INT64 long long int
#endif

#define BOOL INT8
#define TRUE 1
#define FALSE 0

#define PI 3.1415926

#ifndef WIN32
#define stricmp(s1,s2) strcasecmp(s1,s2)
#define strcmpi(s1,s2) strcasecmp(s1,s2)
#endif

#define SQUARE(x)		((x)*(x))
#define CUBE(x)			((x)*(x)*(x))
#define POWER6(x)		SQUARE(CUBE(x))
#define POWER12(x)	SQUARE(POWER6(x))

#endif

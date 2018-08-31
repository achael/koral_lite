/*************************************************************************

                        Mathematica source file

        Copyright 1986 through 1999 by Wolfram Research Inc.


*************************************************************************/

/* C language definitions for use with Mathematica output */

#define Power(x, y)	(pow((ldouble)(x), (ldouble)(y)))
#define Sqrt(x)		(sqrt((ldouble)(x)))
#define Sqrtl(x)        (sqrt((ldouble)(x)))

#define Abs(x)		(fabs((ldouble)(x)))

#define Exp(x)		(exp((ldouble)(x)))
#define Log(x)		(log((ldouble)(x)))

#define Sin(x)		(sin((ldouble)(x)))
#define Cos(x)		(cos((ldouble)(x)))
#define Tan(x)		(tan((ldouble)(x)))

#define ArcSin(x)       (asin((ldouble)(x)))
#define ArcCos(x)       (acos((ldouble)(x)))
#define ArcTan(x)       (atan((ldouble)(x)))

#define Sinh(x)          (sinh((ldouble)(x)))
#define Cosh(x)          (cosh((ldouble)(x)))
#define Tanh(x)          (tanh((ldouble)(x)))

#define Cot(x)          (1./tan((ldouble)(x)))
#define Csc(x)          (1./sin((ldouble)(x)))
#define Sec(x)          (1./cos((ldouble)(x)))

/** Could add definitions for Random(), SeedRandom(), etc. **/



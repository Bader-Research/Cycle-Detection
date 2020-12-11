#include <stdio.h>
#include <stdlib.h>
#include "timing.h"


#ifdef TIME_GTOD
double get_seconds()
{
        struct timeval t;
        struct timezone z;
        gettimeofday(&t,&z);
        return (double)t.tv_sec+((double)t.tv_usec/(double)1000000.0);
}
#endif

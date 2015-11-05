#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sys.h"

#define MA_VERSION "r1"

int main(int argc, char *argv[])
{
	int i;
	sys_init();

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, MA_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, sys_realtime(), sys_cputime());
	return 0;
}

#include "miniasm.h"

int ma_verbose = 3;

void ma_opt_init(ma_opt_t *opt)
{
	opt->min_span = 1000;
	opt->min_match = 100;
	opt->min_dp = 3;
	opt->min_iden = .05;
}

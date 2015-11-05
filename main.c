#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kvec.h"
#include "sys.h"
#include "paf.h"
#include "sdict.h"
#include "miniasm.h"

#define MA_VERSION "r1"

int main(int argc, char *argv[])
{
	ma_opt_t opt;
	int i, c;
	sdict_t *d;
	char *s;
	ma_sub_t *sub;
	ma_hit_t *hit;
	size_t n_hits;

	ma_opt_init(&opt);
	while ((c = getopt(argc, argv, "m:s:d:")) >= 0) {
		if (c == 'm') {
			opt.min_match = strtol(optarg, &s, 10);
			if (*s == ',') opt.min_iden = strtod(s + 1, &s);
		} else if (c == 's') opt.min_span = atoi(optarg);
		else if (c == 'd') opt.min_dp = atoi(optarg);
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: miniasm [options] <in.paf>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  Preselection:\n");
		fprintf(stderr, "    -m INT[,FLOAT]   min match length and fraction [%d,%.2f]\n", opt.min_match, opt.min_iden);
		fprintf(stderr, "    -s INT           min span [%d]\n", opt.min_span);
		fprintf(stderr, "    -d INT           min read depth [%d]\n", opt.min_dp);
		return 1;
	}

	sys_init();
	d = sd_init();

	hit = ma_hit_read(argv[optind], &opt, d, &n_hits);
	sub = ma_hit_sub(opt.min_dp, n_hits, hit);
	n_hits = ma_hit_cut(sub, opt.min_span, n_hits, hit);
	hit = (ma_hit_t*)realloc(hit, n_hits * sizeof(ma_hit_t));

	for (i = 0; i < n_hits; ++i) {
		ma_hit_t *p = &hit[i];
		ma_sub_t *rq = &sub[p->qns>>32], *rt = &sub[p->tn];
		printf("%s:%d-%d\t%d\t%d\t%d\t%c\t%s:%d-%d\t%d\t%d\t%d\n", d->seq[p->qns>>32].name, rq->s + 1, rq->e, rq->e - rq->s, (uint32_t)p->qns, p->qe,
				"+-"[p->rev], d->seq[p->tn].name, rt->s + 1, rt->e, rt->e - rt->s, p->ts, p->te);
	}

	free(sub);
	free(hit);
	sd_destroy(d);

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, MA_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, sys_realtime(), sys_cputime());
	return 0;
}

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

static inline int paf_rec_isflt(const paf_rec_t *r, int min_span, int min_match, float min_frac)
{
	if (r->qe - r->qs < min_span || r->te - r->ts < min_span) return 1;
	if (r->ml < min_match || r->ml < r->bl * min_frac) return 1;
	return 0;
}

int main(int argc, char *argv[])
{
	ma_opt_t opt;
	int i, c;
	paf_file_t *fp;
	sdict_t *d;
	paf_rec_t r;
	char *s;
	ma_hit_v h = {0,0,0};
	ma_reg_t *reg;
	size_t tot = 0;

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
		fprintf(stderr, "  Read pre-selection:\n");
		fprintf(stderr, "    -m INT[,FLOAT]   min match length and fraction [%d,%.2f]\n", opt.min_match, opt.min_iden);
		fprintf(stderr, "    -s INT           min span [%d]\n", opt.min_span);
		fprintf(stderr, "    -d INT           min read depth [%d]\n", opt.min_dp);
		return 1;
	}

	sys_init();
	d = sd_init();

	fp = paf_open(argv[optind]);
	while (paf_read(fp, &r) >= 0) {
		ma_hit_t *p;
		++tot;
		if (paf_rec_isflt(&r, opt.min_span, opt.min_match, opt.min_iden))
			continue;
		kv_pushp(ma_hit_t, h, &p);
		p->qns = (uint64_t)sd_put(d, r.qn, r.ql)<<32 | r.qs;
		p->qe = r.qe;
		p->tn = sd_put(d, r.tn, r.tl);
		p->ts = r.ts, p->te = r.te, p->rev = r.rev;
	}
	paf_close(fp);
	if (ma_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] read %ld hits; stored %ld hits and %d sequences\n", __func__, sys_timestamp(), tot, h.n, d->n_seq);

	ma_hit_sort(h.n, h.a);
	reg = (ma_reg_t*)malloc(d->n_seq * sizeof(ma_reg_t));
	tot = ma_hit_cut(d, opt.min_dp, h.n, h.a, 0, reg);
	if (ma_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] %ld sequences remain after coverage-based cut\n", __func__, sys_timestamp(), tot);

	free(h.a);
	sd_destroy(d);

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, MA_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, sys_realtime(), sys_cputime());
	return 0;
}

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

static void print_subs(const sdict_t *d, const ma_sub_t *sub)
{
	uint32_t i;
	for (i = 0; i < d->n_seq; ++i)
		if (!d->seq[i].del && sub[i].s != sub[i].e)
			printf("%s\t%d\t%d\n", d->seq[i].name, sub[i].s, sub[i].e);
}

static void print_hits(size_t n_hits, const ma_hit_t *hit, const sdict_t *d, const ma_sub_t *sub)
{
	size_t i;
	for (i = 0; i < n_hits; ++i) {
		const ma_hit_t *p = &hit[i];
		const ma_sub_t *rq = &sub[p->qns>>32], *rt = &sub[p->tn];
		printf("%s:%d-%d\t%d\t%d\t%d\t%c\t%s:%d-%d\t%d\t%d\t%d\t%d\t%d\t255\n", d->seq[p->qns>>32].name, rq->s + 1, rq->e, rq->e - rq->s, (uint32_t)p->qns, p->qe,
				"+-"[p->rev], d->seq[p->tn].name, rt->s + 1, rt->e, rt->e - rt->s, p->ts, p->te, p->ml, p->bl);
	}
}

int main(int argc, char *argv[])
{
	ma_opt_t opt;
	int i, c, stage = 100;
	sdict_t *d;
	ma_sub_t *sub = 0;
	ma_hit_t *hit;
	size_t n_hits;
	float cov;
	char *fn_reads = 0, *outfmt = 0;

	ma_opt_init(&opt);
	while ((c = getopt(argc, argv, "m:s:c:S:i:d:g:o:h:I:r:f:e:p:")) >= 0) {
		if (c == 'm') opt.min_match = atoi(optarg);
		else if (c == 'i') opt.min_iden = atof(optarg);
		else if (c == 's') opt.min_span = atoi(optarg);
		else if (c == 'c') opt.min_dp = atoi(optarg);
		else if (c == 'o') opt.min_ovlp = atoi(optarg);
		else if (c == 'S') stage = atoi(optarg);
		else if (c == 'd') opt.bub_dist = atoi(optarg);
		else if (c == 'g') opt.gap_fuzz = atoi(optarg);
		else if (c == 'h') opt.max_hang = atoi(optarg);
		else if (c == 'I') opt.int_frac = atof(optarg);
		else if (c == 'r') opt.ovlp_drop_ratio = atof(optarg);
		else if (c == 'e') opt.max_ext = atoi(optarg);
		else if (c == 'f') fn_reads = optarg;
		else if (c == 'p') outfmt = optarg;
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: miniasm [options] <in.paf>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  Pre-selection:\n");
		fprintf(stderr, "    -m INT      min match length [%d]\n", opt.min_match);
		fprintf(stderr, "    -i FLOAT    min identity [%.2f]\n", opt.min_iden);
		fprintf(stderr, "    -s INT      min span [%d]\n", opt.min_span);
		fprintf(stderr, "    -c INT      min coverage [%d]\n", opt.min_dp);
		fprintf(stderr, "  Overlap:\n");
		fprintf(stderr, "    -o INT      min overlap [%d]\n", opt.min_ovlp);
		fprintf(stderr, "    -h INT      max over hang length [%d]\n", opt.max_hang);
		fprintf(stderr, "    -I FLOAT    min end-to-end match ratio [%.2f]\n", opt.int_frac);
		fprintf(stderr, "  Layout:\n");
		fprintf(stderr, "    -g INT      max gap differences between reads for trans-reduction [%d]\n", opt.gap_fuzz);
		fprintf(stderr, "    -d INT      max distance for bubble popping [%d]\n", opt.bub_dist);
		fprintf(stderr, "    -r FLOAT    overlap drop ratio [%.2f]\n", opt.ovlp_drop_ratio);
		fprintf(stderr, "    -e INT      small unitig threshold [%d]\n", opt.max_ext);
		fprintf(stderr, "    -f FILE     read sequences []\n");
		return 1;
	}

	sys_init();
	d = sd_init();

	hit = ma_hit_read(argv[optind], opt.min_span, opt.min_match, d, &n_hits);

	// first-round filtering
	if (stage >= 2) {
		sub = ma_hit_sub(opt.min_dp, opt.min_iden, 0, n_hits, hit, d->n_seq);
		n_hits = ma_hit_cut(sub, opt.min_span, n_hits, hit);
	}
	if (stage >= 3) n_hits = ma_hit_flt(sub, &opt, n_hits, hit, &cov);

	// second-round filtering
	if (stage >= 4) {
		ma_sub_t *sub2;
		sub2 = ma_hit_sub((int)(cov * .15 + .499) - 1, opt.min_iden, opt.min_span/2, n_hits, hit, d->n_seq);
		n_hits = ma_hit_cut(sub2, opt.min_span, n_hits, hit);
		ma_sub_merge(d->n_seq, sub, sub2);
		free(sub2);
	}
	if (stage >= 5) n_hits = ma_hit_contained(&opt, d, sub, n_hits, hit);
	hit = (ma_hit_t*)realloc(hit, n_hits * sizeof(ma_hit_t));

	// assembly
	if (outfmt != 0 && strcmp(outfmt, "bed") == 0) {
		print_subs(d, sub);
	} else if (outfmt != 0 && strcmp(outfmt, "paf") == 0) {
		print_hits(n_hits, hit, d, sub);
	} if (outfmt == 0 || strcmp(outfmt, "ug") == 0 || strcmp(outfmt, "sg") == 0) {
		asg_t *sg = 0;
		ma_ug_t *ug = 0;

		sg = ma_sg_gen(&opt, d, sub, n_hits, hit);
		asg_arc_del_trans(sg, opt.gap_fuzz);
		asg_pop_bubble(sg, opt.bub_dist);
		asg_cut_short_utg(sg, opt.max_ext, 1);
		asg_arc_del_short(sg, opt.ovlp_drop_ratio);
		asg_pop_bubble(sg, opt.bub_dist);
		for (i = 0; i < 3; ++i)
			asg_cut_short_utg(sg, opt.max_ext, 1);

		if (outfmt == 0 || strcmp(outfmt, "ug") == 0) {
			ug = ma_ug_gen(sg);
			if (fn_reads) ma_ug_seq(ug, d, sub, fn_reads);
			ma_ug_print(ug, d, sub, stdout);
		} else if (strcmp(outfmt, "sg") == 0) {
			ma_sg_print(sg, d, sub, stdout);
		}

		asg_destroy(sg);
		ma_ug_destroy(ug);
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

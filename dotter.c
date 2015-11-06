#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h>
#include "paf.h"
#include "sdict.h"
#include "kvec.h"
#include "eps.h"

typedef struct {
	uint32_t qn, qs, qe;
	uint32_t tn, ts, te;
} dt_hit_t;

typedef struct {
	const char *name;
	uint32_t i;
} srtaux_t;

static inline int mixed_numcompare(const char *_a, const char *_b)
{
	const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
	const unsigned char *pa = a, *pb = b;
	while (*pa && *pb) {
		if (isdigit(*pa) && isdigit(*pb)) {
			while (*pa == '0') ++pa;
			while (*pb == '0') ++pb;
			while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
			if (isdigit(*pa) && isdigit(*pb)) {
				int i = 0;
				while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
				return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
			} else if (isdigit(*pa)) return 1;
			else if (isdigit(*pb)) return -1;
			else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
		} else {
			if (*pa != *pb) return (int)*pa - (int)*pb;
			++pa; ++pb;
		}
	}
	return *pa? 1 : *pb? -1 : 0;
}

#include "ksort.h"
#define srtaux_lt(a, b) (mixed_numcompare((a).name, (b).name) < 0)
KSORT_INIT(dt, srtaux_t, srtaux_lt)

int main(int argc, char *argv[])
{
	int min_span = 1000, min_match = 40, width = 600, height;
	int color[2] = { 0xFF0000, 0x0080FF }, font_size = 11;
	float min_iden = .1;
	paf_file_t *f;
	sdict_t *d[2];
	paf_rec_t r;
	int32_t c, i, j;
	uint64_t *acclen[2], totlen[2];
	srtaux_t *a[2];
	kvec_t(dt_hit_t) h = {0,0,0};
	double sx, sy;

	while ((c = getopt(argc, argv, "m:i:s:w:f:")) >= 0) {
		if (c == 'm') min_match = atoi(optarg);
		else if (c == 'i') min_iden = atof(optarg);
		else if (c == 's') min_span = atoi(optarg);
		else if (c == 'w') width = atoi(optarg);
		else if (c == 'f') font_size = atoi(optarg);
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: minidot [options] <in.paf>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -m INT      min match length [%d]\n", min_match);
		fprintf(stderr, "  -i FLOAT    min identity [%.2f]\n", min_iden);
		fprintf(stderr, "  -s INT      min span [%d]\n", min_span);
		fprintf(stderr, "  -w INT      image width [%d]\n", width);
		fprintf(stderr, "  -f INT      font size [%d]\n", font_size);
		return 1;
	}

	d[0] = sd_init();
	d[1] = sd_init();

	f = paf_open(argv[optind]);
	while (paf_read(f, &r) >= 0) {
		dt_hit_t *s;
		if (r.qe - r.qs < min_span || r.te - r.ts < min_span || r.ml < min_match) continue;
		if (r.ml < r.bl * min_iden) continue;
		kv_pushp(dt_hit_t, h, &s);
		s->qn = sd_put(d[0], r.qn, r.ql), s->qs = r.qs, s->qe = r.qe;
		s->tn = sd_put(d[1], r.tn, r.tl);
		s->ts = r.rev? r.te : r.ts, s->te = r.rev? r.ts : r.te;
	}
	paf_close(f);

	for (i = 0; i < 2; ++i) {
		uint32_t n = d[i]->n_seq;
		uint64_t l = 0;
		a[i] = (srtaux_t*)calloc(n + 1, sizeof(srtaux_t));
		for (j = 0; j < n; ++j)
			a[i][j].name = d[i]->seq[j].name, a[i][j].i = j;
		ks_introsort_dt(d[i]->n_seq, a[i]);
		acclen[i] = (uint64_t*)calloc(n, 8);
		for (j = 0; j < n; ++j)
			acclen[i][a[i][j].i] = l, l += d[i]->seq[a[i][j].i].len;
		totlen[i] = l;
	}
	height = (int)((double)width / totlen[0] * totlen[1] + .499);
	sx = (double)width / totlen[0];
	sy = (double)height / totlen[1];

	eps_header(stdout, width, height, .2);
	eps_font(stdout, "Helvetica-Narrow", font_size);

	// write x label
	eps_gray(stdout, .8);
	for (i = 0; i < d[0]->n_seq; ++i)
		eps_Mstr(stdout, (acclen[0][a[0][i].i] + .5 * d[0]->seq[a[0][i].i].len) * sx, font_size*.5, a[0][i].name);
	eps_stroke(stdout);
	fprintf(stdout, "gsave %g 0 translate 90 rotate\n", font_size*1.25);
	for (i = 0; i < d[1]->n_seq; ++i)
		eps_Mstr(stdout, (acclen[1][a[1][i].i] + .5 * d[1]->seq[a[1][i].i].len) * sx, 0, a[1][i].name);
	fprintf(stdout, "grestore\n");
	eps_stroke(stdout);

	// write grid lines
	eps_linewidth(stdout, .1);
	for (i = 0; i < d[1]->n_seq; ++i)
		eps_linex(stdout, 1, width, i == 0? 1 : acclen[1][a[1][i].i] * sy);
	eps_linex(stdout, 1, width, totlen[1] * sy);
	for (i = 0; i < d[0]->n_seq; ++i)
		eps_liney(stdout, 1, height, i == 0? 1 : acclen[0][a[0][i].i] * sx);
	eps_liney(stdout, 1, height, totlen[0] * sx);
	eps_stroke(stdout);

	// write hits
	eps_linewidth(stdout, .1);
	for (j = 0; j < 2; ++j) {
		eps_color(stdout, color[j]);
		for (i = 0; i < h.n; ++i) {
			dt_hit_t *p = &h.a[i];
			double x0, y0, x1, y1;
			uint64_t xo = acclen[0][p->qn], yo = acclen[1][p->tn];
			if (j == 0 && p->ts > p->te) continue;
			if (j == 1 && p->ts < p->te) continue;
			x0 = (p->qs + xo) * sx, y0 = (p->ts + yo) * sy;
			x1 = (p->qe + xo) * sx, y1 = (p->te + yo) * sy;
			eps_line(stdout, x0, y0, x1, y1);
		}
		eps_stroke(stdout);
	}
	eps_bottom(stdout);

	for (i = 0; i < 2; ++i) {
		free(acclen[i]);
		free(a[i]);
		sd_destroy(d[i]);
	}

	free(h.a);
	return 0;
}

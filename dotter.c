#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include "paf.h"
#include "sdict.h"
#include "kvec.h"
#include "eps.h"

typedef struct {
	uint32_t qn, qs, qe;
	uint32_t tn, ts, te;
	uint32_t ml;
} dt_hit_t;

typedef struct {
	const char *name;
	double tot;
	uint64_t w;
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
#define srtx_lt(a, b) (mixed_numcompare((a).name, (b).name) < 0)
KSORT_INIT(dtx, srtaux_t, srtx_lt)
#define srty_lt(a, b) ((a).tot < (b).tot)
KSORT_INIT(dty, srtaux_t, srty_lt)

int main(int argc, char *argv[])
{
	int min_span = 1000, min_match = 100, width = 600, height, diagonal = 1;
	int color[2] = { 0xFF0000, 0x0080FF }, font_size = 11, no_label = 0;
	float min_iden = .1, lw = 3.0f;
	paf_file_t *f;
	sdict_t *d[2];
	paf_rec_t r;
	int32_t c, i, j;
	uint64_t *acclen[2], totlen[2];
	srtaux_t *a[2];
	kvec_t(dt_hit_t) h = {0,0,0};
	double sx, sy;

	while ((c = getopt(argc, argv, "m:i:s:w:f:Ldt:")) >= 0) {
		if (c == 'm') min_match = atoi(optarg);
		else if (c == 'i') min_iden = atof(optarg);
		else if (c == 's') min_span = atoi(optarg);
		else if (c == 'w') width = atoi(optarg);
		else if (c == 'f') font_size = atoi(optarg);
		else if (c == 'L') no_label = 1;
		else if (c == 'd') diagonal = 0;
		else if (c == 't') lw = atof(optarg);
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: minidot [options] <in.paf>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -m INT      min match length [%d]\n", min_match);
		fprintf(stderr, "  -i FLOAT    min identity [%.2f]\n", min_iden);
		fprintf(stderr, "  -s INT      min span [%d]\n", min_span);
		fprintf(stderr, "  -w INT      image width [%d]\n", width);
		fprintf(stderr, "  -f INT      font size [%d]\n", font_size);
		fprintf(stderr, "  -t FLOAT    line width [%g]\n", lw);
		fprintf(stderr, "  -L          don't print labels\n");
		fprintf(stderr, "  -d          don't try to put hits onto the diagonal\n");
		return 1;
	}

	d[0] = sd_init();
	d[1] = sd_init();

	f = paf_open(argv[optind]);
	if (!f) {
		fprintf(stderr, "[E::%s] could not open PAF file %s\n", __func__, argv[optind]);
		return 1;
	}
	while (paf_read(f, &r) >= 0) {
		dt_hit_t *s;
		if (r.qe - r.qs < min_span || r.te - r.ts < min_span || r.ml < min_match) continue;
		if (r.ml < r.bl * min_iden) continue;
		kv_pushp(dt_hit_t, h, &s);
		s->qn = sd_put(d[1], r.qn, r.ql), s->qs = r.qs, s->qe = r.qe;
		s->tn = sd_put(d[0], r.tn, r.tl);
		s->ts = r.rev? r.te : r.ts, s->te = r.rev? r.ts : r.te;
		s->ml = r.ml;
	}
	paf_close(f);

	for (i = 0; i < 2; ++i) { // 0 for target; 1 for query
		uint32_t n = d[i]->n_seq;
		uint64_t l = 0;
		a[i] = (srtaux_t*)calloc(n + 1, sizeof(srtaux_t));
		if (i == 0 || !diagonal) {
			for (j = 0; j < n; ++j)
				a[i][j].name = d[i]->seq[j].name, a[i][j].i = j;
			ks_introsort_dtx(n, a[i]);
		} else {
			srtaux_t *b = a[i];
			for (j = 0; j < n; ++j)
				b[j].name = d[i]->seq[j].name, b[j].tot = b[j].w = 0, b[j].i = j;
			for (j = 0; j < h.n; ++j) {
				uint64_t w, coor;
				dt_hit_t *p = &h.a[j];
				srtaux_t *q = &b[p->qn];
				coor = acclen[0][p->tn] + (p->ts + p->te) / 2;
				w = (uint64_t)(.01 * p->ml * p->ml + .499);
				q->tot += (double)coor * w;
				q->w += w;
			}
			for (j = 0; j < n; ++j) b[j].tot /= b[j].w;
			ks_introsort_dty(n, b);
		}
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
	eps_gray(stdout, .8);

	if (!no_label) {
		// write x labels
		for (i = 0; i < d[0]->n_seq; ++i)
			eps_Mstr(stdout, (acclen[0][a[0][i].i] + .5 * d[0]->seq[a[0][i].i].len) * sx, font_size*.5, a[0][i].name);
		eps_stroke(stdout);
		fprintf(stdout, "gsave %g 0 translate 90 rotate\n", font_size*1.25);
		// write y labels
		for (i = 0; i < d[1]->n_seq; ++i)
			eps_Mstr(stdout, (acclen[1][a[1][i].i] + .5 * d[1]->seq[a[1][i].i].len) * sx, 0, a[1][i].name);
		fprintf(stdout, "grestore\n");
		eps_stroke(stdout);
	}

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
	eps_linewidth(stdout, lw);
	eps_linecap(stdout, 1);
	for (j = 0; j < 2; ++j) {
		eps_color(stdout, color[j]);
		for (i = 0; i < h.n; ++i) {
			dt_hit_t *p = &h.a[i];
			double x0, y0, x1, y1;
			uint64_t xo = acclen[0][p->tn], yo = acclen[1][p->qn];
			if (j == 0 && p->ts > p->te) continue;
			if (j == 1 && p->ts < p->te) continue;
			x0 = (p->ts + xo) * sx, y0 = (p->qs + yo) * sy;
			x1 = (p->te + xo) * sx, y1 = (p->qe + yo) * sy;
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

#ifndef MINIASM_H
#define MINIASM_H

#include <stdint.h>
#include <sys/types.h>
#include "sdict.h"
#include "asg.h"

extern int ma_verbose;

typedef struct {
	int min_span;
	int min_match;
	int min_dp;
	float min_iden;

	int max_hang;
	int min_ovlp;
	int gap_fuzz;
	float int_frac;
} ma_opt_t;

typedef struct {
	uint64_t qns;
	uint32_t tn:31, cont:1;
	uint32_t ts:31, rev:1;
	uint32_t qe, te;
} ma_hit_t;

typedef struct { size_t n, m; ma_hit_t *a; } ma_hit_v;

typedef struct {
	uint32_t s:31, del:1, e;
} ma_sub_t;

typedef struct {
	sdict_t *d;
	asg_t *g;
} ma_sg_t;

typedef struct {
	uint32_t len:31, circ:1; // len: length of the unitig; circ: circular if non-zero
	uint32_t start, end; // start: starting vertex in the string graph; end: ending vertex
	uint32_t m, n; // number of reads
	uint64_t *a; // list of reads
	char *s; // unitig sequence is not null
} ma_utg_t;

typedef struct { size_t n, m; ma_utg_t *a; } ma_utg_v;

typedef struct {
	ma_utg_v u;
	asg_t *g;
} ma_ug_t;

#ifdef __cplusplus
extern "C" {
#endif

void ma_opt_init(ma_opt_t *opt);
ma_hit_t *ma_hit_read(const char *fn, const ma_opt_t *opt, sdict_t *d, size_t *n);
ma_sub_t *ma_hit_sub(int min_dp, size_t n, const ma_hit_t *a);
size_t ma_hit_cut(const ma_sub_t *reg, int min_span, size_t n, ma_hit_t *a);

#ifdef __cplusplus
}
#endif

#define MA_HT_INT        (-1)
#define MA_HT_QCONT      (-2)
#define MA_HT_TCONT      (-3)
#define MA_HT_SHORT_OVLP (-4)

static inline int ma_hit2arc(const ma_hit_t *h, int ql, int tl, int max_hang, float int_frac, int min_ovlp, asg_arc_t *p)
{
	int32_t tl5, tl3, ext5, ext3, qs = (int32_t)h->qns;
	uint32_t u, v, l; // u: query end; v: target end; l: length from u to v
	u = v = l = UINT32_MAX;
	if (h->rev) tl5 = tl - h->te, tl3 = h->ts; // tl5: 5'-end overhang (on the query strand); tl3: similar
	else tl5 = h->ts, tl3 = tl - h->te;
	ext5 = qs < tl5? qs : tl5;
	ext3 = ql - h->qe < tl3? ql - h->qe : tl3;
	if (ext5 > max_hang || ext3 > max_hang || h->qe - qs < (h->qe - qs + ext5 + ext3) * int_frac)
		return MA_HT_INT;
	if (qs <= tl5 && ql - h->qe <= tl3) return MA_HT_QCONT; // query contained
	else if (qs >= tl5 && ql - h->qe >= tl3) return MA_HT_TCONT; // target contained
	else if (qs > tl5) u = 0, v = !!h->rev, l = qs - tl5;
	else u = 1, v = !h->rev, l = (ql - h->qe) - tl3;
	if (h->qe - qs + ext5 + ext3 < min_ovlp || h->te - h->ts + ext5 + ext3 < min_ovlp) return MA_HT_SHORT_OVLP; // short overlap
	u |= h->qns>>32<<1, v |= h->tn<<1;
	p->ul = (uint64_t)u<<32 | l, p->v = v, p->ol = ql - l, p->del = 0;
	return l;
}

#endif

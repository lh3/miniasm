#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include "sdict.h"
#include "paf.h"
#include "kvec.h"
#include "sys.h"
#include "miniasm.h"

#include "ksort.h"
#define ma_hit_key(a) ((a).qns)
KRADIX_SORT_INIT(hit, ma_hit_t, ma_hit_key, 8)

KSORT_INIT_GENERIC(uint32_t)

void ma_hit_sort(size_t n, ma_hit_t *a)
{
	radix_sort_hit(a, a + n);
}

static inline int paf_rec_isflt(const paf_rec_t *r, int min_span, int min_match, float min_frac)
{
	if (r->qe - r->qs < min_span || r->te - r->ts < min_span) return 1;
	if (r->ml < min_match || r->ml < r->bl * min_frac) return 1;
	return 0;
}

ma_hit_t *ma_hit_read(const char *fn, const ma_opt_t *opt, sdict_t *d, size_t *n)
{
	paf_file_t *fp;
	paf_rec_t r;
	ma_hit_v h = {0,0,0};
	size_t i, tot = 0, tot_len = 0;

	fp = paf_open(fn);
	while (paf_read(fp, &r) >= 0) {
		ma_hit_t *p;
		++tot;
		if (paf_rec_isflt(&r, opt->min_span, opt->min_match, opt->min_iden))
			continue;
		kv_pushp(ma_hit_t, h, &p);
		p->qns = (uint64_t)sd_put(d, r.qn, r.ql)<<32 | r.qs;
		p->qe = r.qe;
		p->tn = sd_put(d, r.tn, r.tl);
		p->ts = r.ts, p->te = r.te, p->rev = r.rev;
	}
	paf_close(fp);
	for (i = 0; i < d->n_seq; ++i)
		tot_len += d->seq[i].len;
	if (ma_verbose >= 3) fprintf(stderr, "[M::%s::%s] read %ld hits; stored %ld hits and %d sequences (%ld bp)\n", __func__, sys_timestamp(), tot, h.n, d->n_seq, tot_len);
	ma_hit_sort(h.n, h.a);
	if (ma_verbose >= 3) fprintf(stderr, "[M::%s::%s] sorted hits\n", __func__, sys_timestamp());
	*n = h.n;
	return h.a;
}

ma_sub_t *ma_hit_sub(int min_dp, size_t n, const ma_hit_t *a)
{
	size_t i, j, last, n_sub = 0, n_remained = 0;
	kvec_t(uint32_t) b = {0,0,0};
	kvec_t(uint64_t) c = {0,0,0};
	ma_sub_t *sub = 0;

	for (i = 0; i < n; ++i) {
		n_sub = n_sub > a[i].qns>>32? n_sub : a[i].qns>>32;
		n_sub = n_sub > a[i].tn? n_sub : a[i].tn;
	}
	++n_sub;
	sub = (ma_sub_t*)calloc(n_sub, sizeof(ma_sub_t));
	for (i = 1, last = 0; i <= n; ++i) {
		if (i == n || a[i].qns>>32 != a[i-1].qns>>32) { // we come to a new query sequence
			size_t start = 0;
			int dp, max_len, max_pos = -1, qid = a[i-1].qns>>32;
			kv_resize(uint32_t, b, i - last);
			b.n = 0;
			for (j = last; j < i; ++j) { // collect all starts and ends
				if (a[j].tn == qid) continue; // skip self match
				kv_push(uint32_t, b, (uint32_t)a[j].qns<<1);
				kv_push(uint32_t, b, a[j].qe<<1|1);
			}
			ks_introsort_uint32_t(b.n, b.a);
			for (j = c.n = 0, dp = 0, max_len = 0; j < b.n; ++j) {
				int old_dp = dp;
				if (b.a[j]&1) --dp;
				else ++dp;
				if (old_dp < min_dp && dp >= min_dp) {
					start = b.a[j]>>1;
				} else if (old_dp >= min_dp && dp < min_dp) {
					int len = (b.a[j]>>1) - start;
					if (max_len < len) max_len = len, max_pos = c.n;
					kv_push(uint64_t, c, (uint64_t)start<<32 | (b.a[j]>>1));
				}
			}
			if (max_pos >= 0) {
				sub[qid].s = (uint32_t)(c.a[max_pos]>>32);
				sub[qid].e = (uint32_t)c.a[max_pos];
				sub[qid].del = 0;
				++n_remained;
			} else sub[qid].del = 1, sub[qid].s = sub[qid].e = 0;
			last = i;
		}
	}
	free(c.a); free(b.a);
	if (ma_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] %ld query sequences remain after sub\n", __func__, sys_timestamp(), n_remained);
	return sub;
}

size_t ma_hit_cut(const ma_sub_t *reg, int min_span, size_t n, ma_hit_t *a)
{
	size_t i, m;
	for (i = m = 0; i < n; ++i) {
		ma_hit_t *p = &a[i];
		const ma_sub_t *rq = &reg[p->qns>>32], *rt = &reg[p->tn];
		int qs, qe, ts, te;
		if (rq->del || rt->del) continue;
		if (p->rev) {
			qs = p->te < rt->e? (uint32_t)p->qns : (uint32_t)p->qns + (p->te - rt->e);
			qe = p->ts > rt->s? p->qe : p->qe - (rt->s - p->ts);
			ts = p->qe < rq->e? p->ts : p->ts + (p->qe - rq->e);
			te = (uint32_t)p->qns > rq->s? p->te : p->te - (rq->s - (uint32_t)p->qns);
		} else {
			qs = p->ts > rt->s? (uint32_t)p->qns : (uint32_t)p->qns + (rt->s - p->ts);
			qe = p->te < rt->e? p->qe : p->qe - (p->te - rt->e);
			ts = (uint32_t)p->qns > rq->s? p->ts : p->ts + (rq->s - (uint32_t)p->qns);
			te = p->qe < rq->e? p->te : p->te - (p->qe - rq->e);
		}
		qs = (qs > rq->s? qs : rq->s) - rq->s;
		qe = (qe < rq->e? qe : rq->e) - rq->s;
		ts = (ts > rt->s? ts : rt->s) - rt->s;
		te = (te < rt->e? te : rt->e) - rt->s;
		if (qe - qs >= min_span && te - ts >= min_span) {
			p->qns = p->qns>>32<<32 | qs, p->qe = qe, p->ts = ts, p->te = te;
			a[m++] = *p;
		}
	}
	if (ma_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] %ld hits remain after cut\n", __func__, sys_timestamp(), m);
	return m;
}

size_t ma_hit_flt(const ma_sub_t *sub, const ma_opt_t *opt, size_t n, ma_hit_t *a, float *cov)
{
	size_t i, m;
	asg_arc_t t;
	uint64_t tot_dp = 0, tot_len = 0;
	int r, max_hang = opt->max_hang + (opt->max_hang>>1), min_ovlp = opt->min_ovlp / 2;
	for (i = m = 0; i < n; ++i) {
		ma_hit_t *h = &a[i];
		const ma_sub_t *sq = &sub[h->qns>>32], *st = &sub[h->tn];
		if (sq->del || st->del) continue;
		r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, .5, min_ovlp, &t);
		if (r >= 0 || r == MA_HT_QCONT || r == MA_HT_TCONT)
			a[m++] = *h, tot_dp += r >= 0? r : r == MA_HT_QCONT? sq->e - sq->s : st->e - st->s;
	}
	for (i = 1; i <= m; ++i)
		if (i == m || a[i].qns>>32 != a[i-1].qns>>32)
			tot_len += sub[a[i-1].qns>>32].e - sub[a[i-1].qns>>32].s;
	*cov = (double)tot_dp / tot_len;
	if (ma_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] %ld hits remain after filtering; crude depth: %.2f\n", __func__, sys_timestamp(), m, *cov);
	return m;
}

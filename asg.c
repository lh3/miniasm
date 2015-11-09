#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "asg.h"
#include "kvec.h"

#include "ksort.h"
#define asg_arc_key(a) ((a).ul)
KRADIX_SORT_INIT(asg, asg_arc_t, asg_arc_key, 8)

asg_t *asg_init(void)
{
	return (asg_t*)calloc(1, sizeof(asg_t));
}

void asg_destroy(asg_t *g)
{
	if (g == 0) return;
	free(g->seq); free(g->idx); free(g->arc); free(g);
}

void asg_arc_sort(asg_t *g)
{
	radix_sort_asg(g->arc, g->arc + g->n_arc);
}

uint64_t *asg_arc_index_core(size_t max_seq, size_t n, const asg_arc_t *a)
{
	size_t i, last;
	uint64_t *idx;
	idx = (uint64_t*)calloc(max_seq * 2, 8);
	for (i = 1, last = 0; i <= n; ++i)
		if (i == n || a[i-1].ul>>32 != a[i].ul>>32)
			idx[a[i-1].ul>>32] = (uint64_t)last<<32 | (i - last), last = i;
	return idx;
}

void asg_arc_index(asg_t *g)
{
	if (g->idx) free(g->idx);
	g->idx = asg_arc_index_core(g->n_seq, g->n_arc, g->arc);
}

void asg_seq_set(asg_t *g, int sid, int len, int del)
{
	if (sid >= g->m_seq) {
		g->m_seq = sid + 1;
		kv_roundup32(g->m_seq);
		g->seq = (asg_seq_t*)realloc(g->seq, g->m_seq * sizeof(asg_seq_t));
	}
	if (sid >= g->n_seq) g->n_seq = sid + 1;
	g->seq[sid].len = len;
	g->seq[sid].del = !!del;
}

// hard remove arcs marked as "del"
void asg_arc_rm(asg_t *g)
{
	uint32_t e, n;
	for (e = n = 0; e < g->n_arc; ++e) {
		uint32_t u = g->arc[e].ul>>32, v = g->arc[e].v;
		if (!g->arc[e].del && !g->seq[u>>1].del && !g->seq[v>>1].del)
			g->arc[n++] = g->arc[e];
	}
	if (n < g->n_arc) { // arc index is out of sync
		if (g->idx) free(g->idx);
		g->idx = 0;
	}
	g->n_arc = n;
}

void asg_cleanup(asg_t *g)
{
	asg_arc_rm(g);
	if (!g->is_srt) {
		asg_arc_sort(g);
		g->is_srt = 1;
	}
	if (g->idx == 0) asg_arc_index(g);
}

// delete short arcs
int asg_arc_del_short(asg_t *g, float drop_ratio)
{
	uint32_t v, n_vtx = g->n_seq * 2, n_short = 0;
	for (v = 0; v < n_vtx; ++v) {
		asg_arc_t *av = asg_arc_a(g, v);
		uint32_t i, thres, nv = asg_arc_n(g, v);
		if (nv < 2) continue;
		thres = (uint32_t)(av[0].ol * drop_ratio + .499);
		for (i = nv - 1; i >= 1 && av[i].ol < thres; --i);
		for (i = i + 1; i < nv; ++i)
			av[i].del = 1, ++n_short;
	}
	if (n_short) {
		asg_cleanup(g);
		asg_symm(g);
	}
	fprintf(stderr, "[M::%s] removed %d short overlaps\n", __func__, n_short);
	return n_short;
}

// delete multi-arcs
int asg_arc_del_multi(asg_t *g)
{
	uint32_t *cnt, n_vtx = g->n_seq * 2, n_multi = 0, v;
	cnt = (uint32_t*)calloc(n_vtx, 4);
	for (v = 0; v < n_vtx; ++v) {
		asg_arc_t *av = asg_arc_a(g, v);
		int32_t i, nv = asg_arc_n(g, v);
		if (nv < 2) continue;
		for (i = nv - 1; i >= 0; --i) ++cnt[av[i].v];
		for (i = nv - 1; i >= 0; --i)
			if (--cnt[av[i].v] != 0)
				av[i].del = 1, ++n_multi;
	}
	free(cnt);
	if (n_multi) asg_cleanup(g);
	fprintf(stderr, "[M::%s] removed %d multi-arcs\n", __func__, n_multi);
	return n_multi;
}

// remove asymmetric arcs: u->v is present, but v'->u' not
int asg_arc_del_asymm(asg_t *g)
{
	uint32_t e, n_asymm = 0;
	for (e = 0; e < g->n_arc; ++e) {
		uint32_t v = g->arc[e].v^1, u = g->arc[e].ul>>32^1;
		uint32_t i, nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		for (i = 0; i < nv; ++i)
			if (av[i].v == u) break;
		if (i == nv) g->arc[e].del = 1, ++n_asymm;
	}
	if (n_asymm) asg_cleanup(g);
	fprintf(stderr, "[M::%s] removed %d asymmetric arcs\n", __func__, n_asymm);
	return n_asymm;
}

void asg_symm(asg_t *g)
{
	asg_arc_del_multi(g);
	asg_arc_del_asymm(g);
	g->is_symm = 1;
}

// transitive reduction; see Myers, 2005
int asg_arc_del_trans(asg_t *g, int fuzz)
{
	uint8_t *mark;
	uint32_t v, n_vtx = g->n_seq * 2, n_reduced = 0;

	mark = (uint8_t*)calloc(n_vtx, 1);
	for (v = 0; v < n_vtx; ++v) {
		uint32_t L, i, nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		if (nv == 0) continue; // no hits
		if (g->seq[v>>1].del) {
			for (i = 0; i < nv; ++i) av[i].del = 1, ++n_reduced;
			continue;
		}
		for (i = 0; i < nv; ++i) mark[av[i].v] = 1;
		L = asg_arc_len(av[nv-1]) + fuzz;
		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].v;
			uint32_t j, nw = asg_arc_n(g, w);
			asg_arc_t *aw = asg_arc_a(g, w);
			if (mark[av[i].v] != 1) continue;
			for (j = 0; j < nw && asg_arc_len(aw[j]) + asg_arc_len(av[i]) <= L; ++j)
				if (mark[aw[j].v]) mark[aw[j].v] = 2;
		}
		#if 0
		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].v;
			uint32_t j, nw = asg_arc_n(g, w);
			asg_arc_t *aw = asg_arc_a(g, w);
			for (j = 0; j < nw && (j == 0 || asg_arc_len(aw[j]) < fuzz); ++j)
				if (mark[aw[j].v]) mark[aw[j].v] = 2;
		}
		#endif
		for (i = 0; i < nv; ++i) {
			if (mark[av[i].v] == 2) av[i].del = 1, ++n_reduced;
			mark[av[i].v] = 0;
		}
	}
	free(mark);
	fprintf(stderr, "[M::%s] transitively reduced %d arcs\n", __func__, n_reduced);
	if (n_reduced) {
		asg_cleanup(g);
		asg_symm(g);
	}
	return n_reduced;
}

/**********************************
 * Filter short potential unitigs *
 **********************************/

#define ASG_ET_MERGEABLE 0
#define ASG_ET_TIP       1
#define ASG_ET_MULTI_OUT 2
#define ASG_ET_MULTI_NEI 3

static inline int asg_is_utg_end(const asg_t *g, uint32_t v, uint64_t *lw)
{
	uint32_t w, nv, nw, nw0, nv0 = asg_arc_n(g, v^1);
	int i, i0 = -1;
	asg_arc_t *aw, *av = asg_arc_a(g, v^1);
	for (i = nv = 0; i < nv0; ++i)
		if (!av[i].del) i0 = i, ++nv;
	if (nv == 0) return ASG_ET_TIP; // tip
	if (nv > 1) return ASG_ET_MULTI_OUT; // multiple outgoing arcs
	if (lw) *lw = av[i0].ul<<32 | av[i0].v;
	w = av[i0].v ^ 1;
	nw0 = asg_arc_n(g, w);
	aw = asg_arc_a(g, w);
	for (i = nw = 0; i < nw0; ++i)
		if (!aw[i].del) ++nw;
	if (nw != 1) return ASG_ET_MULTI_NEI;
	return ASG_ET_MERGEABLE;
}

int asg_extend(const asg_t *g, uint32_t v, int max_ext, asg64_v *a)
{
	int ret;
	uint64_t lw;
	a->n = 0;
	kv_push(uint64_t, *a, v);
	do {
		ret = asg_is_utg_end(g, v^1, &lw);
		if (ret != 0) break;
		kv_push(uint64_t, *a, lw);
		v = (uint32_t)lw;
	} while (--max_ext > 0);
	return ret;
}

int asg_cut_tip(asg_t *g, int max_ext)
{
	asg64_v a = {0,0,0};
	uint32_t n_vtx = g->n_seq * 2, v, i, cnt = 0;
	for (v = 0; v < n_vtx; ++v) {
		if (g->seq[v>>1].del) continue;
		if (asg_is_utg_end(g, v, 0) != ASG_ET_TIP) continue; // not a tip
		if (asg_extend(g, v, max_ext, &a) == ASG_ET_MERGEABLE) continue; // not a short unitig
		for (i = 0; i < a.n; ++i)
			asg_seq_del(g, (uint32_t)a.a[i]>>1);
		++cnt;
	}
	free(a.a);
	if (cnt > 0) asg_cleanup(g);
	fprintf(stderr, "[M::%s] cut %d tips\n", __func__, cnt);
	return cnt;
}

int asg_cut_internal(asg_t *g, int max_ext)
{
	asg64_v a = {0,0,0};
	uint32_t n_vtx = g->n_seq * 2, v, i, cnt = 0;
	for (v = 0; v < n_vtx; ++v) {
		if (g->seq[v>>1].del) continue;
		if (asg_is_utg_end(g, v, 0) != ASG_ET_MULTI_NEI) continue;
		if (asg_extend(g, v, max_ext, &a) != ASG_ET_MULTI_NEI) continue;
		for (i = 0; i < a.n; ++i)
			asg_seq_del(g, (uint32_t)a.a[i]>>1);
		++cnt;
	}
	free(a.a);
	if (cnt > 0) asg_cleanup(g);
	fprintf(stderr, "[M::%s] cut %d internal sequences\n", __func__, cnt);
	return cnt;
}

int asg_cut_biloop(asg_t *g, int max_ext)
{
	asg64_v a = {0,0,0};
	uint32_t n_vtx = g->n_seq * 2, v, i, cnt = 0;
	for (v = 0; v < n_vtx; ++v) {
		uint32_t nv, nw, w, x, ov = 0, ox = 0;
		asg_arc_t *av, *aw;
		if (g->seq[v>>1].del) continue;
		if (asg_is_utg_end(g, v, 0) != ASG_ET_MULTI_NEI) continue;
		if (asg_extend(g, v, max_ext, &a) != ASG_ET_MULTI_OUT) continue;
		x = (uint32_t)a.a[a.n - 1] ^ 1;
		nv = asg_arc_n(g, v ^ 1), av = asg_arc_a(g, v ^ 1);
		for (i = 0; i < nv; ++i)
			if (!av[i].del) w = av[i].v ^ 1;
		nw = asg_arc_n(g, w), aw = asg_arc_a(g, w);
		for (i = 0; i < nw; ++i) { // we are looking for: v->...->x', w->v and w->x
			if (aw[i].del) continue;
			if (aw[i].v == x) ox = aw[i].ol;
			if (aw[i].v == v) ov = aw[i].ol;
		}
		if (ov == 0 && ox == 0) continue;
		if (ov > ox) {
			asg_arc_del(g, w, x, 1);
			asg_arc_del(g, x^1, w^1, 1);
			++cnt;
		}
	}
	free(a.a);
	if (cnt > 0) asg_cleanup(g);
	fprintf(stderr, "[M::%s] cut %d small bi-loops\n", __func__, cnt);
	return cnt;
}

/******************
 * Bubble popping *
 ******************/

typedef struct {
	uint32_t p; // the optimal parent vertex
	uint32_t d; // the shortest distance from the initial vertex
	uint32_t c; // max count of reads
	uint32_t r:31, s:1; // r: the number of remaining incoming arc; s: state
} binfo_t;

typedef struct {
	binfo_t *a;
	kvec_t(uint32_t) S; // set of vertices without parents
	kvec_t(uint32_t) T; // set of tips
	kvec_t(uint32_t) b; // visited vertices
	kvec_t(uint32_t) e; // visited edges/arcs
} buf_t;

// count the number of outgoing arcs, excluding reduced arcs
static inline int count_out(const asg_t *g, uint32_t v)
{
	uint32_t i, n, nv = asg_arc_n(g, v);
	const asg_arc_t *av = asg_arc_a(g, v);
	for (i = n = 0; i < nv; ++i)
		if (!av[i].del) ++n;
	return n;
}

// in a resolved bubble, mark unused vertices and arcs as "reduced"
static void asg_bub_backtrack(asg_t *g, uint32_t v0, buf_t *b)
{
	uint32_t i, v;
	assert(b->S.n == 1);
	for (i = 0; i < b->b.n; ++i)
		g->seq[b->b.a[i]>>1].del = 1;
	for (i = 0; i < b->e.n; ++i) {
		asg_arc_t *a = &g->arc[b->e.a[i]];
		a->del = 1;
		asg_arc_del(g, a->v^1, a->ul>>32^1, 1);
	}
	v = b->S.a[0];
	do {
		uint32_t u = b->a[v].p; // u->v
		g->seq[v>>1].del = 0;
		asg_arc_del(g, u, v, 0);
		asg_arc_del(g, v^1, u^1, 0);
		v = u;
	} while (v != v0);
}

// pop bubbles from vertex v0; the graph MJUST BE symmetric: if u->v present, v'->u' must be present as well
static uint64_t asg_bub_pop1(asg_t *g, uint32_t v0, int max_dist, buf_t *b)
{
	uint32_t i, n_pending = 0;
	uint64_t n_pop = 0;
	if (g->seq[v0>>1].del) return 0; // already deleted
	if ((uint32_t)g->idx[v0] < 2) return 0; // no bubbles
	b->S.n = b->T.n = b->b.n = b->e.n = 0;
	b->a[v0].c = b->a[v0].d = 0;
	kv_push(uint32_t, b->S, v0);
	do {
		uint32_t v = kv_pop(b->S), d = b->a[v].d, c = b->a[v].c;
		uint32_t nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		assert(nv > 0);
		for (i = 0; i < nv; ++i) { // loop through v's neighbors
			uint32_t w = av[i].v, l = (uint32_t)av[i].ul; // u->w with length l
			binfo_t *t = &b->a[w];
			if (w == v0) goto pop_reset;
			if (av[i].del) continue;
			kv_push(uint32_t, b->e, (g->idx[v]>>32) + i);
			if (d + l > max_dist) break; // too far
			if (t->s == 0) { // this vertex has never been visited
				kv_push(uint32_t, b->b, w); // save it for revert
				t->p = v, t->s = 1, t->d = d + l;
				t->r = count_out(g, w^1);
				++n_pending;
			} else { // visited before
				if (c + 1 > t->c || (c + 1 == t->c && d + l > t->d)) t->p = v;
				if (c + 1 > t->c) t->c = c + 1;
				if (d + l < t->d) t->d = d + l; // update dist
			}
			assert(t->r > 0);
			if (--(t->r) == 0) {
				uint32_t x = asg_arc_n(g, w);
				if (x) kv_push(uint32_t, b->S, w);
				else kv_push(uint32_t, b->T, w); // a tip
				--n_pending;
			}
		}
		if (i < nv || b->S.n == 0) goto pop_reset;
	} while (b->S.n > 1 || n_pending);
	asg_bub_backtrack(g, v0, b);
	n_pop = 1 | (uint64_t)b->T.n<<32;
pop_reset:
	for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
		binfo_t *t = &b->a[b->b.a[i]];
		t->s = t->c = t->d = 0;
	}
	return n_pop;
}

// pop bubbles
int asg_pop_bubble(asg_t *g, int max_dist)
{
	uint32_t v, n_vtx = g->n_seq * 2;
	uint64_t n_pop = 0;
	buf_t b;
	if (!g->is_symm) asg_symm(g);
	memset(&b, 0, sizeof(buf_t));
	b.a = (binfo_t*)calloc(n_vtx, sizeof(binfo_t));
	for (v = 0; v < n_vtx; ++v) {
		uint32_t i, n_arc = 0, nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		if (nv < 2 || g->seq[v>>1].del) continue;
		for (i = 0; i < nv; ++i) // asg_bub_pop1() may delete some edges/arcs
			if (!av[i].del) ++n_arc;
		if (n_arc > 1)
			n_pop += asg_bub_pop1(g, v, max_dist, &b);
	}
	free(b.a); free(b.S.a); free(b.T.a); free(b.b.a); free(b.e.a);
	if (n_pop) asg_cleanup(g);
	fprintf(stderr, "[M::%s] popped %d bubbles and trimmed %d tips\n", __func__, (uint32_t)n_pop, (uint32_t)(n_pop>>32));
	return n_pop;
}

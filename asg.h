#ifndef ASG_H
#define ASG_H

#include <stdint.h>
#include <stdlib.h>

typedef struct {
	uint64_t ul;
	uint32_t v;
	uint32_t ol:31, del:1;
} asg_arc_t;

typedef struct {
	uint32_t len:31, del:1;
} asg_seq_t;

typedef struct {
	uint32_t m_arc, n_arc:31, is_srt:1;
	asg_arc_t *arc;
	uint32_t m_seq, n_seq:31, is_symm:1;
	asg_seq_t *seq;
	uint64_t *idx;
} asg_t;

typedef struct { size_t n, m; uint64_t *a; } asg64_v;

#define asg_arc_len(arc) ((uint32_t)(arc).ul)
#define asg_arc_n(g, v) ((uint32_t)(g)->idx[(v)])
#define asg_arc_a(g, v) (&(g)->arc[(g)->idx[(v)]>>32])

asg_t *asg_init(void);
void asg_destroy(asg_t *g);
void asg_seq_set(asg_t *g, int sid, int len, int del);
void asg_symm(asg_t *g);
void asg_cleanup(asg_t *g);

int asg_arc_del_short(asg_t *g, float drop_ratio);
int asg_arc_del_trans(asg_t *g, int fuzz);
int asg_cut_short_utg(asg_t *g, int max_ext, int mask);
int asg_pop_bubble(asg_t *g, int max_dist);

// append an arc
static inline asg_arc_t *asg_arc_pushp(asg_t *g)
{
	if (g->n_arc == g->m_arc) {
		g->m_arc = g->m_arc? g->m_arc<<1 : 16;
		g->arc = (asg_arc_t*)realloc(g->arc, g->m_arc * sizeof(asg_arc_t));
	}
	return &g->arc[g->n_arc++];
}

// set asg_arc_t::del for v->w
static inline void asg_arc_del(asg_t *g, uint32_t v, uint32_t w, int del)
{
	uint32_t i, nv = asg_arc_n(g, v);
	asg_arc_t *av = asg_arc_a(g, v);
	for (i = 0; i < nv; ++i)
		if (av[i].v == w) av[i].del = !!del;
}

#endif

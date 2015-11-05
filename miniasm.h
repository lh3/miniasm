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

#endif

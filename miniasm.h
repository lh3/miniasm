#ifndef MINIASM_H
#define MINIASM_H

#include <stdint.h>
#include <sys/types.h>
#include "sdict.h"

extern int ma_verbose;

typedef struct {
	int min_span;
	int min_match;
	int min_dp;
	float min_iden;
} ma_opt_t;

typedef struct {
	uint64_t qns;
	uint32_t qe, tn, ts:31, rev:1, te;
} ma_hit_t;

typedef struct { size_t n, m; ma_hit_t *a; } ma_hit_v;

typedef struct {
	uint32_t s:31, del:1, e;
} ma_reg_t;

#ifdef __cplusplus
extern "C" {
#endif

void ma_opt_init(ma_opt_t *opt);
size_t ma_hit_cut(const sdict_t *pd, int min_dp, size_t n, ma_hit_t *a, const ma_reg_t *in, ma_reg_t *out);
void ma_hit_sort(size_t n, ma_hit_t *a);

#ifdef __cplusplus
}
#endif

#endif

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include "sdict.h"
#include "paf.h"
#include "kvec.h"
#include "miniasm.h"

#include "ksort.h"
#define ma_hit_key(a) ((a).qns)
KRADIX_SORT_INIT(hit, ma_hit_t, ma_hit_key, 8)

KSORT_INIT_GENERIC(uint32_t)

size_t ma_hit_cut(const sdict_t *pd, int min_dp, size_t n, ma_hit_t *a, ma_reg_t *out)
{
	size_t i, j, last, n_remained = 0;
	kvec_t(uint32_t) b = {0,0,0};
	kvec_t(uint64_t) c = {0,0,0};
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
				out[qid].s = (uint32_t)(c.a[max_pos]>>32);
				out[qid].e = (uint32_t)c.a[max_pos];
				out[qid].del = 0;
				++n_remained;
			} else out[qid].del = 1;
			last = i;
		}
	}
	free(c.a); free(b.a);
	return n_remained;
}

void ma_hit_sort(size_t n, ma_hit_t *a)
{
	radix_sort_hit(a, a + n);
}

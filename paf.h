#ifndef PAF_PAF_H
#define PAF_PAF_H

#include <stdint.h>
#include <sys/types.h>

typedef struct {
	void *fp;
	kstring_t buf;
} paf_file_t;

typedef struct {
	const char *qn, *tn; // these point to the input string; NOT allocated
	uint32_t ql, qs, qe, tl, ts, te;
	uint32_t ml:31, rev:1, bl;
} paf_rec_t;

#ifdef __cplusplus
extern "C" {
#endif

paf_file_t *paf_file_open(const char *fn);
int paf_file_close(paf_file_t *pf);
int paf_file_read(paf_file_t *pf, paf_rec_t *r);

#ifdef __cplusplus
}
#endif

static inline int paf_rec_isflt(const paf_rec_t *r, int min_span, int min_match, float min_frac)
{
	if (r->qe - r->qs < min_span || r->te - r->ts < min_span) return 1;
	if (r->ml < min_match || r->ml < r->bl * min_frac) return 1;
	return 0;
}
#endif

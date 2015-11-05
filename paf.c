#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include "paf.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

paf_file_t *paf_file_open(const char *fn)
{
	kstream_t *ks;
	gzFile fp;
	paf_file_t *pf;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	pf = (paf_file_t*)calloc(1, sizeof(paf_file_t));
	pf->fp = ks;
	return pf;
}

int paf_file_close(paf_file_t *pf)
{
	kstream_t *ks;
	if (pf == 0) return 0;
	free(pf->buf.s);
	ks = (kstream_t*)pf->fp;
	gzclose(ks->f);
	ks_destroy(ks);
	free(pf);
	return 0;
}

int paf_parse(int l, char *s, paf_rec_t *pr) // s must be NULL terminated
{ // on return: <0 for failure; 0 for success; >0 for filtered
	char *q, *r;
	int i, t;
	for (i = t = 0, q = s; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) pr->qn = q;
		else if (t == 1) pr->ql = strtol(q, &r, 10);
		else if (t == 2) pr->qs = strtol(q, &r, 10);
		else if (t == 3) pr->qe = strtol(q, &r, 10);
		else if (t == 4) pr->rev = (*q == '-');
		else if (t == 5) pr->tn = q;
		else if (t == 6) pr->tl = strtol(q, &r, 10);
		else if (t == 7) pr->ts = strtol(q, &r, 10);
		else if (t == 8) pr->te = strtol(q, &r, 10);
		else if (t == 9) pr->ml = strtol(q, &r, 10);
		else if (t == 10) pr->bl = strtol(q, &r, 10);
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 10) return -1;
	return 0;
}

int paf_file_read(paf_file_t *pf, paf_rec_t *r)
{
	int ret, dret;
file_read_more:
	ret = ks_getuntil((kstream_t*)pf->fp, KS_SEP_LINE, &pf->buf, &dret);
	if (ret < 0) return ret;
	ret = paf_parse(pf->buf.l, pf->buf.s, r);
	if (ret < 0) goto file_read_more;
	return ret;
}

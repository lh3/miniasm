#include <string.h>
#include "sdict.h"

#include "khash.h"
KHASH_MAP_INIT_STR(str, uint32_t)
typedef khash_t(str) shash_t;

sdict_t *sd_init(void)
{
	sdict_t *d;
	d = (sdict_t*)calloc(1, sizeof(sdict_t));
	d->h = kh_init(str);
	return d;
}

void sd_destroy(sdict_t *d)
{
	uint32_t i;
	if (d == 0) return;
	if (d->h) kh_destroy(str, (shash_t*)d->h);
	for (i = 0; i < d->n_seq; ++i)
		free(d->seq[i].name);
	free(d->seq);
	free(d);
}

int32_t sd_put(sdict_t *d, const char *name, uint32_t len)
{
	shash_t *h = (shash_t*)d->h;
	khint_t k;
	int absent;
	k = kh_put(str, h, name, &absent);
	if (absent) {
		sd_seq_t *s;
		if (d->n_seq == d->m_seq) {
			d->m_seq = d->m_seq? d->m_seq<<1 : 16;
			d->seq = (sd_seq_t*)realloc(d->seq, d->m_seq * sizeof(sd_seq_t));
		}
		s = &d->seq[d->n_seq];
		s->len = len, s->aux = 0, s->del = 0;
		kh_key(h, k) = s->name = strdup(name);
		kh_val(h, k) = d->n_seq++;
	} // TODO: test if len is the same;
	return kh_val(h, k);
}

int32_t *sd_squeeze(sdict_t *d)
{
	int32_t *map, i, j;
	if (d->h) {
		kh_destroy(str, (shash_t*)d->h);
		d->h = 0;
	}
	map = (int32_t*)calloc(d->n_seq, 4);
	for (i = j = 0; i < d->n_seq; ++i) {
		if (d->seq[i].del) {
			free(d->seq[i].name);
			map[i] = -1;
		} else d->seq[j] = d->seq[i], map[i] = j++;
	}
	d->n_seq = j;
	return map;
}

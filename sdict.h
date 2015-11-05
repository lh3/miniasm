#ifndef SDICT_H
#define SDICT_H

#include <stdint.h>

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

typedef struct {
	char *name;
	uint32_t len, aux;
} sd_seq_t;

typedef struct {
	uint32_t n_seq, m_seq;
	sd_seq_t *seq;
	void *h;
} sdict_t;

#ifdef __cplusplus
extern "C" {
#endif

sdict_t *sd_init(void);
void sd_destroy(sdict_t *d);
int32_t sd_put(sdict_t *d, const char *name, uint32_t len);

#ifdef __cplusplus
}
#endif

#endif

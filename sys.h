#ifndef HL_SYS_H
#define HL_SYS_H

#ifdef __cplusplus
extern "C" {
#endif

double sys_cputime();
double sys_realtime();
void sys_init();
const char *sys_timestamp();

#ifdef __cplusplus
}
#endif

#endif

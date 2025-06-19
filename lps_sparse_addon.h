#ifndef _LPS_SPARSE_ADDON_
#define _LPS_SPARSE_ADDON_
#include "lps.h"
int parse1_sparse(const char *begin, const char *end, struct core *cores, uint64_t offset);
void init_sparse_lps_offset(struct lps *lps_ptr, const char *str, int len, uint64_t offset);
#endif

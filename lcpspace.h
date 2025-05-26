#ifndef __LCPSPACE_H__
#define __LCPSPACE_H__

#include "struct_def.h"
#include "utils.h"
#include "tpool.h"
#include <stdio.h>
#include <string.h>

/**
 * @brief Constructs a LCP space DBG and prints in rGFA format.
 *
 * @param seqs       A pointer to the `ref_seq` structure containing the reference
 *                   sequences and their processed LCP cores.
 * @param out        A file pointer to the output file where the formatted segments
 *                   and links will be written.
 */
void lcp_space_graph_construct_from_genome(struct ref_seq *seqs, FILE *out);

#endif

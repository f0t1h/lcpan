#include "lps_sparse_addon.h"

void init_sparse_lps_offset(struct lps *lps_ptr, const char *str, int len, uint64_t offset) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores = (struct core *)malloc((len/CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse1_sparse(str, str+len, lps_ptr->cores, offset);
}


template<class FUNC>
void init_core1_sparse(struct core *cr, const char *begin, uint64_t distance, uint64_t start_index, uint64_t end_index, FUNC cond, int *idx){
    init_core1(cr, begin, distance, start_index, end_index);
    if(cond(cr)){
        ++(*idx);
    }
}

bool core_ok(struct core *c) { return 10000000 > c->label; }

int parse1_sparse(const char *begin, const char *end, struct core *cores, uint64_t offset) {

    const char *it1 = begin;
    const char *it2 = end;
    int core_index = 0;
    int last_invalid_char_index = -1;

    // find lcp cores
    for (; it1 + 2 < end; it1++) {

        // skip invalid character
        if (alphabet[(unsigned char)*it1] == -1) {
            last_invalid_char_index = it1 - begin;
            continue;
        }

        if (alphabet[(unsigned char)*it1] == alphabet[(unsigned char)*(it1+1)]) {
            continue;
        }

        // check for RINT core
        if (alphabet[(unsigned char)*(it1+1)] == alphabet[(unsigned char)*(it1+2)]) {

            // count middle characters
            uint32_t middle_count = 1;
            const char *temp = it1 + 2;
            while (temp < end && alphabet[(unsigned char)*(temp-1)] == alphabet[(unsigned char)*temp]) {
                temp++;
                middle_count++;
            }
            if (temp != end) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                    init_core1_sparse(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset, core_ok, &core_index);
                }

                // create RINT core
                it2 = it1 + 2 + middle_count;
                init_core1_sparse(&(cores[core_index]), it1, it2-it1, it1-begin+offset, it2-begin+offset, core_ok, &core_index);

                continue;
            }
        }

        if (alphabet[(unsigned char)*it1] > alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] < alphabet[(unsigned char)*(it1+2)]) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                init_core1_sparse(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset, core_ok, &core_index);
            }

            // create LMIN core
            it2 = it1 + 3;
            init_core1_sparse(&(cores[core_index]), it1, it2-it1, it1-begin+offset, it2-begin+offset, core_ok, &core_index);

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (it1+3 < end &&
            alphabet[(unsigned char)*it1] < alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] > alphabet[(unsigned char)*(it1+2)] &&
            alphabet[(unsigned char)*(it1-1)] <= alphabet[(unsigned char)*(it1)] &&
            alphabet[(unsigned char)*(it1+2)] >= alphabet[(unsigned char)*(it1+3)]) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                init_core1_sparse(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset, core_ok, &core_index);
            }

            // create LMAX core
            it2 = it1 + 3;
            init_core1_sparse(&(cores[core_index]), it1, it2-it1, it1-begin+offset, it2-begin+offset, core_ok, &core_index);
            continue;
        }
    }

    return core_index;
}


/*
int parse1_sparse(const char *begin, const char *end, struct core *cores, uint64_t offset) {

    const char *it1 = begin;
    const char *it2 = end;
    int core_index = 0;
    int last_invalid_char_index = -1;

    // find lcp cores
    for (; it1 + 2 < end; it1++) {

        // skip invalid character
        if (alphabet[(unsigned char)*it1] == -1) {
            last_invalid_char_index = it1 - begin;
            continue;
        }

        if (alphabet[(unsigned char)*it1] == alphabet[(unsigned char)*(it1+1)]) {
            continue;
        }

        // check for RINT core

        if (begin == it1) {
            continue;
        }
#ifndef USE_LMINS
#define USE_LMINS 0
#endif
#ifndef KEEP_LMIN_CORES
#define KEEP_LMIN_CORES 1
#endif
#if USE_LMINS
        if (alphabet[(unsigned char)*it1] > alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] < alphabet[(unsigned char)*(it1+2)]) {



            // create LMIN core
            it2 = it1 + 3;
#if KEEP_LMIN_CORES
            init_core1(&(cores[core_index]), it1, it2-it1, it1-begin+offset, it2-begin+offset);
            if(core_index == 0 || cores[core_index].bit_rep != cores[core_index-1].bit_rep){ // Skip identical adjacent cores
                core_index++;
            }
#endif
            continue;
        }
#endif
        // check for LMAX
        if (it1+3 < end &&
            alphabet[(unsigned char)*it1] < alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] > alphabet[(unsigned char)*(it1+2)] &&
            alphabet[(unsigned char)*(it1-1)] <= alphabet[(unsigned char)*(it1)] 
            && alphabet[(unsigned char)*(it1+2)] >= alphabet[(unsigned char)*(it1+3)]
            ) {
            // create LMAX core
            it2 = it1 + 3;
            init_core1(&(cores[core_index]), it1, it2-it1, it1-begin+offset, it2-begin+offset);
            //if(core_index == 0 || cores[core_index].bit_rep != cores[core_index-1].bit_rep){ // Skip identical adjacent cores
            core_index++;
            //}

            continue;
        }
    }

    return core_index;
}
*/


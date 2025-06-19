#include <algorithm>
#include <unordered_map>
#include <array>
#include <vector>
#include <map>
#include "fa_parser.h"

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "lbdg.h"

#define MIN_REQUIRED_CORE_COUNT 2
#define LCMER_SIZE 8
#define TIME_CHECKPOINT_INIT(NAME)\
	clock_t NAME##_start = clock(), NAME##_diff;\
	int NAME##_msec
	
#define TIME_CHECKPOINT(NAME, fmt) do{\
		NAME##_diff = clock() - NAME##_start;\
		NAME##_msec = NAME##_diff * 1000 / CLOCKS_PER_SEC;\
		fprintf(stderr, fmt, NAME##_msec/1000, NAME##_msec%1000);\
		NAME##_start = clock();\
	}while(0)
struct lcmer{
    union {
        std::array<uint32_t, LCMER_SIZE> ar;
    }data;

    bool operator <(const lcmer &other)const {
        for(int i = 0; i < LCMER_SIZE; ++i){
            if(data.ar[i] != other.data.ar[i]){
                if(data.ar[i] < other.data.ar[i]){
                    return true;
                }
                else{
                    return false;
                }
            }
        }
        return false;
    }
    void print_name(FILE *out) const{
        fprintf(out, "%0x", data.ar[0]);
        for(int z = 1; z < LCMER_SIZE; ++z){
            fprintf(out,"%0x", data.ar[z]);
        }
    }
};
void lcmer_set0(lcmer &l, uint32_t new_core){
    l.data.ar[LCMER_SIZE-1] = new_core;
}
lcmer shift_lcmer(lcmer l){
    for(int i = 1; i < 8; ++i){
        l.data.ar[i-1] = l.data.ar[i];
    }
    return l;
}
lcmer update_lcmer(lcmer l, uint32_t new_core){
    for(int i = 1; i < 8; ++i){
        l.data.ar[i-1] = l.data.ar[i];
    }
    l.data.ar[LCMER_SIZE-1] = new_core;
    return l;
}

void lspag_print_ref_seq(struct opt_arg *args, FILE *out) {
    
    TIME_CHECKPOINT_INIT(LSPACETIME);
    std::unordered_map<uint64_t, uint64_t> core_counter;
    printf("[INFO] Processing reference...\n");

    fprintf(out, "H\tVN:Z:1.1\n");
/*

    for (int i=0; i<seqs->size; i++) {
        const struct chr chrom = seqs->chrs[i];
        for (int j=0; j<chrom.cores_size; j++) core_counter[chrom.cores[j].id]++;
    }
    


    for (int i=0; i<seqs->size; i++) {
        struct chr &chrom = seqs->chrs[i];
        int right = 0;
        int left;
        for(left = 0; left < chrom.cores_size; ++left){
            while(right < chrom.cores_size && core_counter[chrom.cores[right].id] < MIN_REQUIRED_CORE_COUNT){
                ++right;
            }
            if(right == chrom.cores_size) { break;}
            if(right > left){
                chrom.cores[left] = chrom.cores[right];
                ++right;
            }
        }
        chrom.cores_size = left;
    }
*/

    std::map<lcmer, std::vector<uint32_t>> dbg;
   
    gzFile fp;
    kseq_t *seq;


    fp = gzopen(args->fasta_path, "r");
    seq = kseq_init(fp);
    size_t idx = 0;
    lcmer l;
    while(kseq_read(seq) >= 0){
        struct chr chrom = {seq->name.s, idx++, (int) seq->seq.l, seq->seq.s, 0, 0, 0};
        lbdg_process_chrom(seq->seq.s, seq->seq.l, args->lcp_level, &chrom);

        struct simple_core *cores = chrom.cores;
        if (chrom.cores_size < LCMER_SIZE + 2){
            continue;
        }
        for(int j = 1; j < LCMER_SIZE + 1; ++j){
            l.data.ar[j-1] = cores[j].id;   
        }
        for (int j=LCMER_SIZE+1; j<chrom.cores_size - 1; j++) {
            dbg[l].push_back(cores[j].id);
            l = update_lcmer(l, cores[j].id);
        }
        dbg[l];
        free(chrom.cores);
    }
    kseq_destroy(seq);

    gzclose(fp);
    /*
    for (int i=0; i<seqs->size; i++) {
        struct chr &chrom = seqs->chrs[i];
		lcmer l;
        struct simple_core *cores = chrom.cores;
        if (chrom.cores_size < LCMER_SIZE + 2){
            continue;
        }
        for(int j = 1; j < LCMER_SIZE + 1; ++j){
            l.data.ar[j-1] = cores[j].id;   
        }
        for (int j=LCMER_SIZE+1; j<chrom.cores_size - 1; j++) {
            dbg[l].push_back(cores[j].id);
            l = update_lcmer(l, cores[j].id);
        }
        dbg[l];
	}
    */
    TIME_CHECKPOINT(LSPACETIME, "Built DBG %d sec %d ms\n");
    for(const auto &p : dbg){
        fprintf(out,"S\t");
        p.first.print_name(out);
        fprintf(out,"\t*\tLN:i:1000\n");
    }
    for(const auto &p :dbg){
        lcmer ln = shift_lcmer(p.first);
        for(uint32_t next : p.second){
            fprintf(out,"L\t");
            p.first.print_name(out);
            fprintf(out,"\t+\t");
            lcmer_set0(ln, next);
            ln.print_name(out);
            fprintf(out,"\t+\t100M\n");
        }
    }
    TIME_CHECKPOINT(LSPACETIME, "Printed DBG %d sec %d ms\n");

}


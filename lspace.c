#include <algorithm>
#include <unordered_map>
#include <array>
#include <vector>
#include <map>

#include "lbdg.h"

#define MIN_REQUIRED_CORE_COUNT 2
#define LCMER_SIZE 8
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
lcmer update_lcmer(lcmer l, uint32_t new_core){
    for(int i = 1; i < 8; ++i){
        l.data.ar[i-1] = l.data.ar[i];
    }
    l.data.ar[LCMER_SIZE-1] = new_core;
    return l;
}

void lspag_print_ref_seq(struct ref_seq *seqs, FILE *out) {
    

    std::unordered_map<uint64_t, uint64_t> core_counter;
    printf("[INFO] Processing reference...\n");

    fprintf(out, "H\tVN:Z:1.1\n");


    for (int i=0; i<seqs->size; i++) {
        const struct chr chrom = seqs->chrs[i];
        for (int j=0; j<chrom.cores_size; j++) core_counter[chrom.cores[j].id]++;
    }
    

/*
    for (int i=0; i<seqs->size; i++) {
        struct chr &chrom = seqs->chrs[i];
        auto end = std::remove_if( chrom.cores, chrom.cores + chrom.cores_size,
                [&core_counter] (uint64_t e) -> bool{
                    return core_counter[e] < MIN_REQUIRED_CORE_COUNT;
                    }
                );
        chrom.cores_size = end - chrom.cores;
    }
*/

    std::map<lcmer, std::vector<uint32_t>> dbg;

    for (int i=0; i<seqs->size; i++) {
        struct chr &chrom = seqs->chrs[i];
		lcmer l;
        struct simple_core *cores = chrom.cores;
        if (chrom.cores_size < 10){
            continue;
        }
        for(int j = 0; j < 8; ++j){
            l.data.ar[j] = cores[j].id;   
        }
        for (int j=8; j<chrom.cores_size; j++) {
            dbg[l].push_back(cores[j].id);
            l = update_lcmer(l, cores[j].id);
        }
	}
    
    for(const auto &p : dbg){
        fprintf(out,"S\t");
        p.first.print_name(out);
        fprintf(out,"\t*\tLN:i:1000\n");
    }
    for(const auto &p :dbg){
        for(uint32_t next : p.second){
            fprintf(out,"L\t");
            p.first.print_name(out);
            fprintf(out,"\t+\t");
            update_lcmer(p.first, next).print_name(out);
            fprintf(out,"\t+\t100M\n");
        }
    }

}


#include <string.h>
#include <zlib.h>

#include <algorithm>
#include <array>
#include <map>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "fa_parser.h"
#include "include/blockingconcurrentqueue.h"
#include "include/gtl/phmap.hpp"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "lbdg.h"

#define MIN_REQUIRED_CORE_COUNT 2
#define LCMER_SIZE 24
#define TIME_CHECKPOINT_INIT(NAME)               \
    clock_t NAME##_start = clock(), NAME##_diff; \
    int NAME##_msec

#define TIME_CHECKPOINT(NAME, fmt)                                    \
    do {                                                              \
        NAME##_diff = clock() - NAME##_start;                         \
        NAME##_msec = NAME##_diff * 1000 / CLOCKS_PER_SEC;            \
        fprintf(stderr, fmt, NAME##_msec / 1000, NAME##_msec % 1000); \
        NAME##_start = clock();                                       \
    } while (0)
struct lcmer {
    union {
        std::array<uint32_t, LCMER_SIZE> ar;
    } data;

    bool operator==(const lcmer &other) const {
        for (int i = 0; i < LCMER_SIZE; ++i) {
            if (data.ar[i] != other.data.ar[i]) return false;
        }
        return true;
    }
    bool operator<(const lcmer &other) const {
        for (int i = 0; i < LCMER_SIZE; ++i) {
            if (data.ar[i] != other.data.ar[i]) {
                if (data.ar[i] < other.data.ar[i]) {
                    return true;
                } else {
                    return false;
                }
            }
        }
        return false;
    }
    void print_name(FILE *out) const {
        fprintf(out, "%0x", data.ar[0]);
        for (int z = 1; z < LCMER_SIZE; ++z) {
            fprintf(out, "%0x", data.ar[z]);
        }
    }
};
void lcmer_set0(lcmer &l, uint32_t new_core) {
    l.data.ar[LCMER_SIZE - 1] = new_core;
}
lcmer shift_lcmer(lcmer l) {
    for (int i = 1; i < LCMER_SIZE; ++i) {
        l.data.ar[i - 1] = l.data.ar[i];
    }
    return l;
}
lcmer update_lcmer(lcmer l, uint32_t new_core) {
    for (int i = 1; i < LCMER_SIZE; ++i) {
        l.data.ar[i - 1] = l.data.ar[i];
    }
    l.data.ar[LCMER_SIZE - 1] = new_core;
    return l;
}

template <>
struct std::hash<lcmer> {
    std::size_t operator()(const lcmer &l) const noexcept {
        return MurmurHash3_32(&l.data.ar, LCMER_SIZE * sizeof(uint32_t));
    }
};

void lspag_print_ref_seq(struct opt_arg *args, FILE *out) {
    TIME_CHECKPOINT_INIT(LSPACETIME);

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
       while(right < chrom.cores_size && core_counter[chrom.cores[right].id] <
       MIN_REQUIRED_CORE_COUNT){
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

    //    std::map<lcmer, std::vector<uint32_t>> dbg;
    using mutex_type = std::mutex;
    using next_lcmer_set =
        gtl::parallel_flat_hash_set<uint32_t,
                                    gtl::priv::hash_default_hash<uint32_t>,
                                    gtl::priv::hash_default_eq<uint32_t>,
                                    std::allocator<uint32_t>, 4, mutex_type>;
    using dbg_map = gtl::parallel_flat_hash_map<
        lcmer, next_lcmer_set, std::hash<lcmer>,
        gtl::priv::hash_default_eq<lcmer>,
        std::allocator<std::pair<const lcmer, next_lcmer_set>>, 4, mutex_type>;

    dbg_map dbg;
    gzFile fp;
    kseq_t *seq;

    fp = gzopen(args->fasta_path, "r");
    seq = kseq_init(fp);
    size_t idx = 0;

    auto process_seq = [&](char *n, char *s, int len) {
        lcmer l;
        struct chr chrom = COLITERAL(chr){n, idx++, len, s, 0, 0, 0};
        lbdg_process_chrom(s, len, args->lcp_level, &chrom);
       // printf("%d\n", chrom.cores_size);
        struct simple_core *cores = chrom.cores;
        if (chrom.cores_size < LCMER_SIZE + 2) {
            return;
        }
        for (int j = 1; j <= LCMER_SIZE; ++j) {
            l.data.ar[j - 1] = cores[j].id;
        }
        for (int j = LCMER_SIZE + 1; j < chrom.cores_size - 1; j++) {
            dbg[l].insert(cores[j].id);
            l = update_lcmer(l, cores[j].id);
        }
        dbg[l];
        free(chrom.cores);
        free(n);
        free(s);
    };
    using produced_data_type = std::tuple<char *, char *, int>;
    moodycamel::BlockingConcurrentQueue<produced_data_type> bcq;
    moodycamel::ProducerToken pt(bcq);
    bool done;
    std::thread fqproducer([&]() {
        while (kseq_read(seq) >= 0) {
            char *n = (char *)malloc(seq->name.l + 1);
            char *s = (char *)malloc(seq->seq.l + 1);
            std::strncpy(n, seq->name.s, seq->name.l + 1);
            std::strncpy(s, seq->seq.s, seq->seq.l + 1);
            bcq.enqueue(pt, std::make_tuple(n, s, seq->seq.l));
        }
        done = true;
    });

    std::vector<std::thread> consumers;

    for (int i = 0; i < std::max(1, args->thread_number-1); ++i) {
        consumers.emplace_back([&]() {
            moodycamel::ConsumerToken ct(bcq);
            std::vector<produced_data_type> queried;
            while (!done || bcq.size_approx() > 0) {
                std::tuple<char *, char *, int> item;

                bcq.wait_dequeue_bulk(ct, std::back_inserter(queried), 10);
                for(const auto &item : queried)
                    process_seq(std::get<0>(item),
                            std::get<1>(item),
                            std::get<2>(item));
                queried.clear();
            }
        });
    }
    fqproducer.join();
    for (std::thread &t : consumers) {
        t.join();
    }

    /*
       while(kseq_read(seq) >= 0){
       struct chr chrom = COLITERAL(chr){seq->name.s, idx++, (int) seq->seq.l,
       seq->seq.s, 0, 0, 0}; lbdg_process_chrom(seq->seq.s, seq->seq.l,
       args->lcp_level, &chrom);

       struct simple_core *cores = chrom.cores;
       if (chrom.cores_size < LCMER_SIZE + 2){
       continue;
       }
       for(int j = 1; j <= LCMER_SIZE; ++j){
       l.data.ar[j-1] = cores[j].id;
       }
       for (int j=LCMER_SIZE+1; j<chrom.cores_size-1; j++) {
       dbg[l].insert(cores[j].id);
       l = update_lcmer(l, cores[j].id);
       }
       dbg[l];
       free(chrom.cores);
       }
       */
    kseq_destroy(seq);

    gzclose(fp);

    TIME_CHECKPOINT(LSPACETIME, "Built DBG %d sec %d ms\n");
    for (const auto &p : dbg) {
        fprintf(out, "S\t");
        p.first.print_name(out);
        fprintf(out, "\t*\tLN:i:1000\n");
    }
    for (const auto &p : dbg) {
        lcmer ln = shift_lcmer(p.first);
        for (uint32_t next : p.second) {
            fprintf(out, "L\t");
            p.first.print_name(out);
            fprintf(out, "\t+\t");
            lcmer_set0(ln, next);
            ln.print_name(out);
            fprintf(out, "\t+\t100M\n");
        }
    }
    TIME_CHECKPOINT(LSPACETIME, "Printed DBG %d sec %d ms\n");
}


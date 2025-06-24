
#include <string.h>
#include <zlib.h>

#include <array>
#include <thread>
#include <tuple>
#include <utility>
#include <chrono>
#include <atomic>
#include "include/concurrentqueue.h"

#include "include/gtl/phmap.hpp"
#include "include/gtl/vector.hpp"
#include <stdio.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#include "fa_parser.h"

#ifndef LCMER_SIZE
#define LCMER_SIZE 32
#endif

#include <chrono>
#include <ctime>
#include <cstdio>


#define INIT_TIMERS(name) \
    std::chrono::high_resolution_clock::time_point name##_chrono_start = std::chrono::high_resolution_clock::now(); \
    std::chrono::high_resolution_clock::time_point name##_last_checkpoint = name##_chrono_start; \
    clock_t name##_clock_start = std::clock(); \
    clock_t name##_last_clock = name##_clock_start;


#define TIME_CHECKPOINT(name, stream, fmt, ...) \
    do { \
        auto name##_now = std::chrono::high_resolution_clock::now(); \
        clock_t name##_clock_now = std::clock(); \
        \
        long long name##_total_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(name##_now - name##_chrono_start).count(); \
        long long name##_interval_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(name##_now - name##_last_checkpoint).count(); \
        \
        double name##_total_cpu = 1000.0 * (name##_clock_now - name##_clock_start) / CLOCKS_PER_SEC; \
        double name##_interval_cpu = 1000.0 * (name##_clock_now - name##_last_clock) / CLOCKS_PER_SEC; \
        \
        double name##_total_util = (name##_total_elapsed > 0) ? (name##_total_cpu / name##_total_elapsed) * 100.0 : 0.0; \
        double name##_interval_util = (name##_interval_elapsed > 0) ? (name##_interval_cpu / name##_interval_elapsed) * 100.0 : 0.0; \
        \
        fprintf(stream, "[Checkpoint: %s] ", #name); \
        fprintf(stream, fmt, ##__VA_ARGS__); \
        fprintf(stream, "\n"); \
        fprintf(stream, "  Total:   Wall = %lld ms, CPU = %.2f ms, Utilization = %.2f%%\n", name##_total_elapsed, name##_total_cpu, name##_total_util); \
        fprintf(stream, "  Current: Wall = %lld ms, CPU = %.2f ms, Utilization = %.2f%%\n\n", name##_interval_elapsed, name##_interval_cpu, name##_interval_util); \
        \
        name##_last_checkpoint = name##_now; \
        name##_last_clock = name##_clock_now; \
    } while (0)

static std::atomic<unsigned long long> thread_counter;

unsigned long long thread_id() {
    thread_local unsigned long long tid = ++thread_counter;
    return tid;
}

template<int N>
struct lcmer {
    std::array<uint32_t, N> data;

    lcmer () {}

    bool operator==(const lcmer &other) const {
        for (int i = 0; i < N; ++i) {
            if (data[i] != other.data[i]) return false;
        }
        return true;
    }
    bool operator<(const lcmer &other) const {
        for (int i = 0; i < N; ++i) {
            if (data[i] != other.data[i]) {
                if (data[i] < other.data[i]) {
                    return true;
                } else {
                    return false;
                }
            }
        }
        return false;
    }
    void print_name(FILE *out) const {
        fprintf(out, "%0x", data[0]);
        for (int z = 1; z < N; ++z) {
            fprintf(out, "%0x", data[z]);
        }
    }


};
template<int N>
void lcmer_set0(lcmer<N> &l, uint32_t new_core) {
    l.data[N - 1] = new_core;
}
template<int N>
lcmer<N> shift_lcmer(lcmer<N> l) {
    for (int i = 1; i < N; ++i) {
        l.data[i - 1] = l.data[i];
    }
    return l;
}
template<int N>
lcmer<N> update_lcmer(lcmer<N> l, uint32_t new_core) {
    for (int i = 1; i < N; ++i) {
        l.data[i - 1] = l.data[i];
    }
    l.data[N - 1] = new_core;
    return l;
}


template <int N>
struct std::hash<lcmer<N>> {
    std::size_t operator()(const lcmer<N> &l) const noexcept {
        return MurmurHash3_32(&l.data, N * sizeof(uint32_t));
    }
};

template<int N>
struct lcmerp_eq{
    inline bool operator() (const lcmer<N> *a, const lcmer<N> *b)const {
        return *a==*b;
    }
};

template<int N>
struct lcmerp_hash{
    inline bool operator() (const lcmer<N> *l)const{
        return MurmurHash3_32(&l->data, N * sizeof(uint32_t));
    }

};

template<int N>
struct lcmerp_heq{
    using Hash = lcmerp_hash<N>;
    using Eq   = lcmerp_eq<N>;
};

struct lcrun{
    std::vector<uint32_t> data;
    template<int N>
    lcrun(const lcmer<N> &lc) : data(std::begin(lc.data), std::end(lc.data)){}
    void push_back(uint32_t d) {
        data.push_back(d);
    }
};
template<>
struct std::hash<lcrun> {
    std::size_t operator()(const lcrun &l) const noexcept {
        return MurmurHash3_32(l.data.data(), l.data.size() * sizeof(uint32_t));
    }
};
void lspag_print_ref_seq(struct opt_arg *args, FILE *out) {
    INIT_TIMERS(LSPACETIME);

    fprintf(stderr,"[INFO] Processing reference...\n");

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
    using idx_and_nls = std::pair<size_t, next_lcmer_set>;
    using dbg_map = gtl::parallel_flat_hash_map<
        lcmer<LCMER_SIZE>, idx_and_nls, std::hash<lcmer<LCMER_SIZE>>,
        gtl::priv::hash_default_eq<lcmer<LCMER_SIZE>>,
        std::allocator<std::pair<const lcmer<LCMER_SIZE>, idx_and_nls>>, 4, mutex_type>;

    using core_pos_map_t = gtl::parallel_flat_hash_map<
        uint64_t,
        std::tuple<char *, int, int, int, int, struct simple_core *>,
        std::hash<uint64_t>,
        gtl::priv::hash_default_eq<uint64_t>,
        std::allocator<std::pair<const uint64_t, std::tuple<char*, int, int, int, int, struct simple_core *>>>,
        4,
        mutex_type>;
    core_pos_map_t cpm;

    dbg_map dbg;
    gzFile fp;
    kseq_t *seq;

    fp = gzopen(args->fasta_path, "r");
    seq = kseq_init(fp);
    size_t idx = 0;

    TIME_CHECKPOINT(LSPACETIME, stderr, "Starting Assembly");
    auto process_seq = [&](char *n, char *s, int len, int tid, uint64_t lcmer_index) {
        lcmer<LCMER_SIZE> l;
        struct chr chrom = COLITERAL(chr){n, idx++, len, s, 0, 0, 0};
        lbdg_process_chrom(s, len, args->lcp_level, &chrom);
       // printf("%d\n", chrom.cores_size);
        struct simple_core *cores = chrom.cores;
        bool str_used = false;
        if (chrom.cores_size < LCMER_SIZE + 2) {
            return lcmer_index;
        }
        uint64_t cpair = cores[0].id;
        for (int j = 1; j <= LCMER_SIZE; ++j) {
            l.data[j - 1] = cores[j].id;
            cpair = (cpair << 32) | (unsigned) cores[j].id;
            str_used |= cpm.try_emplace(cpair, s, cores[j].start, cores[j].end, j, chrom.cores_size, cores).second;

        }
        for (int j = LCMER_SIZE + 1; j < chrom.cores_size - 1; j++) {
            dbg.lazy_emplace_l(
                    l,
                    [&](dbg_map::value_type &v){
                        v.second.second.insert(cores[j].id);
                    },
                    [&](const dbg_map::constructor &ctor){
                        ctor(l, std::pair(lcmer_index | tid, next_lcmer_set{}));
                        lcmer_index+=128;
                    }
                    );
            cpair = (cpair << 32) | (unsigned) cores[j].id;
            str_used |= cpm.try_emplace(cpair, s, cores[j].start, cores[j].end, j, chrom.cores_size, cores).second;
//            dbg[l].insert(cores[j].id);
            l = update_lcmer(l, cores[j].id);
        }

        dbg.lazy_emplace_l(
                l,
                [](dbg_map::value_type &){},
                [&](const dbg_map::constructor &ctor){
                    ctor(l, std::pair(lcmer_index | tid, next_lcmer_set{}));
                    lcmer_index += 128;
                }
                );
//        dbg[l];
//        free(chrom.cores);
        free(n);
        if(!str_used){
            free(s);
        }
        return lcmer_index;
    };
    using produced_data_type = std::tuple<char *, char *, int>;
    moodycamel::ConcurrentQueue<produced_data_type> bcq;
    moodycamel::ProducerToken pt(bcq);
    std::atomic<bool> done;
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
            unsigned int tid = thread_id();
            fprintf(stderr,"t %d\n", tid);

            uint64_t idx = 0;
            moodycamel::ConsumerToken ct(bcq);
            gtl::vector<produced_data_type> queried;
            while (!done || bcq.size_approx() > 0) {
                bcq.try_dequeue_bulk(ct, std::back_inserter(queried), 32);

                for(const auto &item : queried)
                    idx = process_seq(std::get<0>(item),
                            std::get<1>(item),
                            std::get<2>(item),
                            tid, idx);
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
       l.data[j-1] = cores[j].id;
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
    
    
    using lcmer_idx_map_t = gtl::parallel_flat_hash_map<
        size_t,const lcmer<LCMER_SIZE> *>;
    
//    std::unordered_map<const lcmer<LCMER_SIZE> *, size_t, lcmerp_hash<LCMER_SIZE>, lcmerp_eq<LCMER_SIZE>> lim;
    lcmer_idx_map_t lim;
    idx=0;
    TIME_CHECKPOINT(LSPACETIME, stderr, "Built DBG: LCMER #=%lu", dbg.size() );
    for (const auto &p : dbg) {
        fprintf(out, "S\t%lu\t*\tLN:i:%d\n", p.second.first, LCMER_SIZE);
        lim[p.second.first] = &p.first;
 //       lim[&p.first] = idx++;
 //       fprintf(out, "S\t");
//        p.first.print_name(out);
 //       fprintf(out, "\t*\tLN:i:%d\n",LCMER_SIZE);
    }
    for (const auto &p : dbg) {
        lcmer ln = shift_lcmer(p.first);
  //      size_t lindex = lim.at(&p.first);
        for (uint32_t next : p.second.second) {
            lcmer_set0(ln, next);
//            size_t rindex = lim.at(&ln);
              fprintf(out,"L\t%lu\t+\t%lu\t+\t%dM\n",p.second.first,dbg[ln].first,LCMER_SIZE-1);
  //          fprintf(out, "L\t");
 //           p.first.print_name(out);
//            fprintf(out, "\t+\t");
//            ln.print_name(out);
//            fprintf(out, "\t+\t%dM\n",LCMER_SIZE-1);
        }
    }
    fclose(out);

    char cmd[1024];
    snprintf(cmd, 1024, "gfatools asm -u %s", args->gfa_path); 
    FILE *gfatools_p = popen(cmd, "r");
    size_t N = 0;
    char *buffer = NULL;
    
    using unitig_map_t = gtl::parallel_flat_hash_map<std::string, lcrun>;
    unitig_map_t um;
    
    using unitig_adj_t = gtl::parallel_flat_hash_map<std::string, gtl::vector<std::string>>;
    unitig_adj_t ua;

    using lcmer_2_unitig_pos_t = gtl::parallel_flat_hash_map<
        const lcmer<LCMER_SIZE> *, std::pair<std::string, size_t>,
        lcmerp_hash<LCMER_SIZE>,
        lcmerp_eq<LCMER_SIZE>,
        std::allocator<std::pair<const lcmer<LCMER_SIZE> *,std::pair<std::string, size_t>>>,
        4>;
    lcmer_2_unitig_pos_t  find_unitigs;

    out = fopen("proc.gfa","w");
    while(getline(&buffer, &N, gfatools_p)!=-1){
        fprintf(out, "%s", buffer);
        char *token = strtok(buffer, "\t");
        switch(*token){
        case 'S':
            break;
        case 'A':{
            token = strtok(NULL,"\t");
            strtok(NULL,"\t");
            std::string uid{token};
            strtok(NULL,"\t");
            token = strtok(NULL,"\t");
            const lcmer<LCMER_SIZE> *l = lim.at(atoi(token));
            um.lazy_emplace_l(
                    uid,
                    [&](unitig_map_t::value_type &v){
//                        find_unitigs.try_emplace(l, uid, v.second.data.size());
                        v.second.push_back(l->data[LCMER_SIZE-1]);
                    },
                    [&](const unitig_map_t::constructor &ctor){
                        ctor(uid, *l);
                    }
            );
            }
            break;
        case 'L':{//TODO Implement orientation logic
            token = strtok(NULL,"\t");
            std::string left {token};
            token = strtok(NULL,"\t");
            token = strtok(NULL,"\t");
            std::string right {token};
            ua[left].push_back(right);
            break;
            }
        default:
            fprintf(stderr, "Unknown type %c\n", *token);
        }
    }

    uint64_t bcount = 0;
    TIME_CHECKPOINT(LSPACETIME, stderr, "Simplified DBG");
    for(const auto &p : um){
        printf(">%s\t%lu\n", p.first.c_str(), p.second.data.size());
        uint64_t core = ((uint64_t)p.second.data[0] << 32) | p.second.data[1];

        auto &core_loc = cpm.at(core);

        int start = std::get<1>(core_loc);
        int end = std::get<2>(core_loc);
        char *sq = std::get<0>(core_loc);
        for(size_t i = 2; i < p.second.data.size(); ++i){
            uint64_t pcore = core;
            core = (core << 32) | p.second.data[i];
            auto &core_loc = cpm.at(core);
            if(sq==std::get<0>(core_loc)){
                start = std::min(start, std::get<1>(core_loc));
                end   = std::max(end,   std::get<2>(core_loc));
            }
            else{
                printf("%.*s", end-start, sq+start);
                bcount+=(end-start);
                start = std::get<1>(core_loc);
                end   = std::get<2>(core_loc);
                //free(sq);
                sq    = std::get<0>(core_loc);
                auto *core_l = std::get<5>(cpm.at(pcore));
                int pp = std::get<3>(cpm.at(pcore));
                if(core_l[pp+1].end < core_l[pp].end){
                    continue;
                }
                int ovlp = core_l[pp].end - core_l[pp+1].start;
                if (ovlp > 0){
                    start += ovlp;
                }

//                uint64_t trio[3] = {core_l[pp-1].id,core_l[pp].id,core_l[pp+1].id};
//                fprintf(stderr, "%lu\t%lu\t%lu\t%lu\n", core, pcore, (trio[0]<<32)|trio[1], (trio[1]<<32)|trio[2]);

            }
        }
        printf("%.*s\n", end-start, sq+start);
        bcount+=end-start;
    }
    //TODO Using longer lcmer to resolve ambiguity
    //TODO iterate reads, replace the cores with sequences using the idx
    fclose(out);
    fclose(gfatools_p);
    free(buffer);
    TIME_CHECKPOINT(LSPACETIME, stderr, "Printed Assembly Sequences %lu unitigs and %llu bases", um.size(), bcount);
}


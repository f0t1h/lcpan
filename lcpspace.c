#include "lcpspace.h"
#include "ringbuffer.h"
#include "fa_parser.h"
#include "m-bptree.h"
#include "m-array.h"
#include "m-list.h"
#include "m-dict.h"
#include <time.h>
uint64_t hash_combine(uint64_t seed, uint64_t val){
    seed ^= val + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}
uint32_t murmurhash (const char *key, uint32_t len, uint32_t seed) {
  uint32_t c1 = 0xcc9e2d51;
  uint32_t c2 = 0x1b873593;
  uint32_t r1 = 15;
  uint32_t r2 = 13;
  uint32_t m = 5;
  uint32_t n = 0xe6546b64;
  uint32_t h = 0;
  uint32_t k = 0;
  uint8_t *d = (uint8_t *) key; // 32 bit extract from `key'
  const uint32_t *chunks = NULL;
  const uint8_t *tail = NULL; // tail - last 8 bytes
  int i = 0;
  int l = len / 4; // chunk length

  h = seed;

  chunks = (const uint32_t *) (d + l * 4); // body
  tail = (const uint8_t *) (d + l * 4); // last 8 byte chunk of `key'

  // for each 4 byte chunk of `key'
  for (i = -l; i != 0; ++i) {
    // next 4 byte chunk of `key'
  #if MURMURHASH_HAS_HTOLE32
    k = htole32(chunks[i]);
  #else
    k = chunks[i];
  #endif

    // encode next 4 byte chunk of `key'
    k *= c1;
    k = (k << r1) | (k >> (32 - r1));
    k *= c2;

    // append to hash
    h ^= k;
    h = (h << r2) | (h >> (32 - r2));
    h = h * m + n;
  }

  k = 0;

  // remainder
  switch (len & 3) { // `len % 4'
    case 3: k ^= (tail[2] << 16);
    __attribute__ ((fallthrough));
    case 2: k ^= (tail[1] << 8);
    __attribute__ ((fallthrough));
    case 1:
      k ^= tail[0];
      k *= c1;
      k = (k << r1) | (k >> (32 - r1));
      k *= c2;
      h ^= k;
  }

  h ^= len;

  h ^= (h >> 16);
  h *= 0x85ebca6b;
  h ^= (h >> 13);
  h *= 0xc2b2ae35;
  h ^= (h >> 16);

  return h;
}
uint64_t hash_64_fnv1a(const void* key, const uint64_t len) {
    
    const char* data = (char*)key;
    uint64_t hash = 0xcbf29ce484222325;
    uint64_t prime = 0x100000001b3;
    
    for(uint64_t i = 0; i < len; ++i) {
        uint8_t value = data[i];
        hash = hash ^ value;
        hash *= prime;
    }
    
    return hash;

} //hash_64_fnv1a

uint64_t FNVHashInt(uint64_t val){
    return hash_64_fnv1a((void *)&val, 8);
}

uint64_t MurmurHash3_int(uint64_t val){
    return  murmurhash((void *) &val, 8, 42);
}

#define HASH_COUNT 2
#define LCMER_SIZE 18
//DEFINE_HASHED_RING_BUFFER(hashrb, HASH_COUNT, LCMER_SIZE, FNVHashInt, MurmurHash3_int, uint64_t) 

struct core_loc {
	int start;
	int end;
	uint64_t tid;
};

struct lcmer_t{
	uint64_t data[LCMER_SIZE];
};

typedef struct lcmer2_t{
	int64_t *from;
	int64_t *to;
}lcmer2_t;

uint64_t lcmer_hash(struct lcmer_t *l){
	uint64_t ret = 0;
	for(int i = 0; i < LCMER_SIZE; ++i){
		ret = hash_combine(ret, FNVHashInt(l->data[i]));
	}
	return ret;	
}
//ARRAY_DEF(lcmer, uint64_t)
ARRAY_DEF(array_core_locs, struct core_loc, M_POD_OPLIST)
#define M_OPL_array_core_locs  ARRAY_OPLIST(array_core_locs, M_POD_OPLIST)
//BPTREE_DEF2(dbg_adj, 4, uint64_t, M_BASIC_OPLIST, array_core_locs_t, M_OPL_array_core_locs)
BPTREE_DEF(dbg_adj, 4, uint64_t, M_BASIC_OPLIST)
#define M_OPL_dbg_adj  BPTREE_OPLIST2(dbg_adj, M_BASIC_OPLIST, M_OPL_array_core_locs)

//BPTREE_DEF2(debruin, 4, hashrb_t, M_OPEXTEND(M_POD_OPLIST, HASH(hashrb_rehash2), EQUAL(hashrb_eq), CMP(hashrb_cmp)), dbg_adj_t, M_OPL_dbg_adj)

//BPTREE_DEF2(debruin, 4, struct lcmer_t, M_POD_OPLIST, dbg_adj_t, M_OPL_dbg_adj)

uint64_t range_hash(lcmer2_t data){
	uint64_t hashval = 0;
	for(int64_t *it = data.from; it != data.to; ++it){
		hashval = hash_combine(hashval, *it);
	}
	return hashval;
}
int range_cmp(lcmer2_t a, lcmer2_t b){
	int64_t *itb = b.from;
	for(int64_t *ita = a.from; ita != a.to; ++ita, ++itb){
		if(*ita!=*itb){
			return (int64_t)(*ita) - *itb;
		}
	}
	return 0;	
}
bool range_equal(lcmer2_t a, lcmer2_t b){
	return range_cmp(a,b) == 0;
}
#define M_RANGE_OPLIST M_OPEXTEND(M_POD_OPLIST, HASH(range_hash), EQUAL(range_equal), CMP(range_cmp))
DICT_DEF2(debruin, struct lcmer2_t, M_RANGE_OPLIST, dbg_adj_t, M_OPL_dbg_adj)
DICT_DEF2(lcmer_indexer, struct lcmer2_t, M_RANGE_OPLIST, uint64_t, M_BASIC_OPLIST)
#define M_OPL_lcmer2_t() M_POD_OPLIST


ARRAY_DEF(array_core_id, uint64_t)


#define SHIFT_COUNT (64/LCMER_SIZE)
uint64_t hash_move_right(uint64_t h, uint64_t n){
	return h ^ n;	
}
uint64_t hash_move_left(uint64_t h, uint64_t n){
	return h ^ n;
}


// Maybe a bloom filter 
#define MIN_REQUIRED_CORE_COUNT 2
BPTREE_DEF2(counter, 4, uint64_t, size_t)

BPTREE_DEF2(lcmer_counter, 32, struct lcmer_t, M_POD_OPLIST, size_t, M_BASIC_OPLIST)


void lcmer_shift_left(struct lcmer_t *l){
	for(int i = 0; i < LCMER_SIZE -1; ++i){
		l->data[i]=l->data[i+1];
	}
}
/*
void dfs_step(debruin_t *r, dbg_adj_t *adj, lcmer_counter_t *c, struct lcmer_t *lc, FILE *out, uint64_t prev_hash){

	size_t *visited = lcmer_counter_safe_get(*c, *lc);
	if(*visited == 0){
		++(*visited);
		dbg_adj_it_t ait;
		lcmer_shift_left(lc);
		for(dbg_adj_it(ait, *adj); !dbg_adj_end_p(ait); dbg_adj_next(ait)){
			uint64_t new_core = *dbg_adj_ref(ait)->key_ptr;
			lc->data[LCMER_SIZE-1] = new_core;
			struct lcmer_t new_lc;
			memcpy(&new_lc, &lc, sizeof(struct lcmer_t)); 
			uint64_t hv = lcmer_hash(&new_lc);
//			fprintf(out, "L\t%lu\t+\t%lu\t+\t*\n", prev_hash, hv);
			dfs_step(r, debruin_safe_get(*r, new_lc), c, &new_lc, out, hv);
		}
	}
}

void dfs(debruin_t *r, lcmer_counter_t *c, FILE *out){
	
	debruin_it_t it;
	for(debruin_it(it, *r); !debruin_end_p(it); debruin_next(it)){
		struct lcmer_t *lc = debruin_ref(it)->key_ptr;
		struct lcmer_t new_lc;
		memcpy(&new_lc, lc, sizeof(struct lcmer_t));
		uint64_t hv = lcmer_hash(&new_lc);
		dfs_step(r, debruin_ref(it)->value_ptr, c, &new_lc, out, hv);
	}

}
*/

#define TIME_CHECKPOINT_INIT(NAME)\
	clock_t NAME##_start = clock(), NAME##_diff;\
	int NAME##_msec;
	
#define TIME_CHECKPOINT(NAME, fmt) do{\
		NAME##_diff = clock() - NAME##_start;\
		NAME##_msec = NAME##_diff * 1000 / CLOCKS_PER_SEC;\
		fprintf(stderr, fmt, NAME##_msec/1000, NAME##_msec%1000);\
		NAME##_start = clock();\
	}while(0)
void lcp_space_graph_construct_from_genome(struct ref_seq *seqs, FILE *out) {

	TIME_CHECKPOINT_INIT(LCPTIME)

	counter_t counter_dict;
	counter_init(counter_dict);

	debruin_t left_dbg;
	debruin_init(left_dbg);

	debruin_t right_dbg;
	debruin_init(right_dbg);
		
	lcmer_indexer_t lci_dic;
	lcmer_indexer_init(lci_dic);



	for(int ch_i = 0; ch_i < seqs->size; ++ch_i){
		int core_i;
		for(core_i = 0; core_i < seqs->chrs[ch_i].cores_size; ++core_i){
			size_t *count = counter_safe_get( counter_dict, seqs->chrs[ch_i].cores[core_i].id);
			++(*count);
		}
	}

	TIME_CHECKPOINT(LCPTIME, "Core counting took %d seconds %d milliseconds\n");	

	//Remove singles
	for(int ch_i = 0; ch_i < seqs->size; ++ch_i){
		fprintf(stderr, "%d size before selection: %lu\n", ch_i, seqs->chrs[ch_i].cores_size); 
		int core_j = 0;
		int core_i;
		for(core_i = 0; core_j < seqs->chrs[ch_i].cores_size; ++core_i, ++core_j){
			size_t *count;
			do{
				count = counter_safe_get( counter_dict, seqs->chrs[ch_i].cores[core_j].id);
				++core_j;
			}while(core_j < seqs->chrs[ch_i].cores_size && *count < MIN_REQUIRED_CORE_COUNT);
			if(core_j == seqs->chrs[ch_i].cores_size){
				break;
			}
			seqs->chrs[ch_i].cores[core_i] = seqs->chrs[ch_i].cores[core_j];
		}
		seqs->chrs[ch_i].cores_size = core_i;
		fprintf(stderr,"%d size after selection: %lu\n", ch_i, seqs->chrs[ch_i].cores_size); 
	}

	TIME_CHECKPOINT(LCPTIME,"Core filtering took %d seconds %d milliseconds\n");	


	uint64_t **core_ids=malloc(sizeof(uint64_t *) * seqs->size);

	for(int ch_i = 0; ch_i < seqs->size; ++ch_i){
		core_ids[ch_i] = malloc(sizeof(uint64_t) * seqs->chrs[ch_i].cores_size);
		for(int core_i = 0; core_i < seqs->chrs[ch_i].cores_size; ++core_i){
			core_ids[ch_i][core_i] =  seqs->chrs[ch_i].cores[core_i].id;
		}
	}
	uint64_t lcmer_idx = 1;
	for(int ch_i = 0; ch_i < seqs->size; ++ch_i){
		uint64_t prev = 0;
		for(int core_i = 0; core_i < seqs->chrs[ch_i].cores_size - LCMER_SIZE; ++core_i){
			uint64_t core_in = core_ids[ch_i][ core_i + LCMER_SIZE];
			lcmer2_t lcmer = {core_ids[ch_i] + core_i, core_ids[ch_i] + core_i+LCMER_SIZE};
			dbg_adj_t *rad = debruin_safe_get(right_dbg, lcmer);

			int start1 = seqs->chrs[ch_i].cores[core_i].start;
			int start2 = seqs->chrs[ch_i].cores[core_i+LCMER_SIZE-1].start;

			int end1   = seqs->chrs[ch_i].cores[core_i].end;
			int end2   = seqs->chrs[ch_i].cores[core_i+LCMER_SIZE-1].end;
			
			uint64_t *idx = lcmer_indexer_safe_get(lci_dic, lcmer);
			if(*idx == 0){
				*idx = lcmer_idx;
				lcmer_idx++;
				fprintf(out, "S\t");
				fprintf(out, "C%lu", *idx);
				//for(int i = 0; i < LCMER_SIZE; ++i){
				//	fprintf(out, "-%lu", lcl.data[i]);
				//}
				fprintf(out, "\t*\tRC:i:%d\tLN:i:%u\n", ch_i, LCMER_SIZE * (end2 - start2));
			}
			if(rad==NULL){
				fprintf(stderr, "Cannot get at %d:%d-%d\n", ch_i, core_i, core_i + LCMER_SIZE);
			}
			else{
				dbg_adj_push(*rad, core_in);
				if(dbg_adj_size(*rad) == 1){// && array_core_locs_size(*dbg_adj_safe_get(*rad, core_in))== 1){

				}
				if(core_i != 0){
					fprintf(out, "L\tC%lu\t+\tC%lu\t+\t%dM\tID:Z:%d\n", prev, *idx, start2 - start1,ch_i);
				}

			}
			prev = *idx;
		}

	}
	/*
	for(int ch_i = 0; ch_i < seqs->size; ++ch_i){
		int core_i = 0;
		if(0) core_i = seqs->chrs[ch_i].cores_size;
		int fill = 0;
		struct lcmer_t lcl = {0};
		
		//while(core_i < seqs->chrs[ch_i].cores_size && fill < LCMER_SIZE){
		if(0) { // reversed
			while(core_i >= 0 && fill < LCMER_SIZE){
				size_t *count = counter_safe_get( counter_dict, seqs->chrs[ch_i].cores[core_i].id);
				if(*count >= MIN_REQUIRED_CORE_COUNT){
					lcl.data[fill] = seqs->chrs[ch_i].cores[core_i].id;
					++fill;
	//				lcmer_push_back(lcl, seqs->chrs[ch_i].cores[core_i].id);

					//hashrb_put( &rb, seqs->chrs[ch_i].cores[core_i].id);
				}
				--core_i;
			}
		}
		else{
			while(core_i < seqs->chrs[ch_i].cores_size && fill < LCMER_SIZE){
				size_t *count = counter_safe_get( counter_dict, seqs->chrs[ch_i].cores[core_i].id);
				if(*count >= MIN_REQUIRED_CORE_COUNT){
					lcl.data[fill] = seqs->chrs[ch_i].cores[core_i].id;
					++fill;
	//				lcmer_push_back(lcl, seqs->chrs[ch_i].cores[core_i].id);

					//hashrb_put( &rb, seqs->chrs[ch_i].cores[core_i].id);
				}
				++core_i;
			}

		}

		uint64_t prev = 0;
		if(0){//reversed
			for(int i = 0; i < LCMER_SIZE/2;++i){
				int tmp = lcl.data[i];
				lcl.data[i] = lcl.data[LCMER_SIZE-1-i];
				lcl.data[LCMER_SIZE-1-i] = tmp;
			}
		      core_i = 0;
		}
		for(; core_i < seqs->chrs[ch_i].cores_size; ++core_i){
			size_t *count = counter_safe_get( counter_dict, seqs->chrs[ch_i].cores[core_i].id);
			if(*count >= MIN_REQUIRED_CORE_COUNT){

				uint64_t core_out = 0;
				uint64_t core_in = seqs->chrs[ch_i].cores[core_i].id;
				struct simple_core sc = seqs->chrs[ch_i].cores[core_i];
				struct core_loc cl = {sc.start, sc.end, ch_i};
				dbg_adj_t *rad = debruin_safe_get(right_dbg, lcl);
				uint64_t hv = lcmer_hash(&lcl);
				if(rad==NULL){
					fprintf(stderr, "Cannot get %lu\n", hv); 
					continue;
				}
//				array_core_locs_push_back(*dbg_adj_safe_get(*rad, core_in), cl);
				dbg_adj_push(*rad, core_in);
				if(dbg_adj_size(*rad) == 1){// && array_core_locs_size(*dbg_adj_safe_get(*rad, core_in))== 1){
					fprintf(out, "S\t");
					fprintf(out, "%lu", hv);
					//for(int i = 0; i < LCMER_SIZE; ++i){
					//	fprintf(out, "-%lu", lcl.data[i]);
					//}
					fprintf(out, "\t*\tRC:i:%d\tLN:i:%lu\n", ch_i, LCMER_SIZE * (sc.end - sc.start));
				}
				if(prev!=0){
					fprintf(out, "L\t%lu\t+\t%lu\t+\t*\tID:Z:%d\n", prev, hv, ch_i);
				}

				for(int i = 0; i < LCMER_SIZE -1; ++i){
					lcl.data[i]=lcl.data[i+1];
				}
				lcl.data[LCMER_SIZE-1] = core_in;
				prev = hv;
			}
		}
	}
*/
	TIME_CHECKPOINT(LCPTIME,"Graph building %d seconds %d milliseconds\n");	

/*
	start = clock();
	lcmer_counter_t cd;
	lcmer_counter_init(cd);
	dfs(&right_dbg, &cd, out);

	diff = clock() - start;
	msec = diff * 1000 / CLOCKS_PER_SEC;
	printf("DFS took %d seconds %d milliseconds\n", msec/1000, msec%1000);	
*/	
	for(int ch_i = 0; ch_i < seqs->size; ++ch_i){
		free(core_ids[ch_i]);
	}
	free(core_ids);
	counter_clear(counter_dict);
	debruin_clear(right_dbg);
	debruin_clear(left_dbg);
}

#include "lcpspace.h"
#include "ringbuffer.h"
#include "fa_parser.h"
#include "m-bptree.h"
#include "m-array.h"
#include "m-list.h"

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
#define BUFFER_SIZE 14
//DEFINE_HASHED_RING_BUFFER(hashrb, HASH_COUNT, BUFFER_SIZE, FNVHashInt, MurmurHash3_int, uint64_t) 

struct core_loc {
	int start;
	int end;
	uint64_t tid;
};

struct lcmer_t{
	uint64_t data[BUFFER_SIZE];
};

uint64_t lcmer_hash(struct lcmer_t *l){
	uint64_t ret = 0;
	for(int i = 0; i < BUFFER_SIZE; ++i){
		ret = hash_combine(ret, FNVHashInt(l->data[i]));
	}
	return ret;	
}
//ARRAY_DEF(lcmer, uint64_t)
ARRAY_DEF(array_core_locs, struct core_loc, M_POD_OPLIST)
#define M_OPL_array_core_locs  ARRAY_OPLIST(array_core_locs, M_POD_OPLIST)
BPTREE_DEF2(dbg_adj, 4, uint64_t, M_BASIC_OPLIST, array_core_locs_t, M_OPL_array_core_locs)
#define M_OPL_dbg_adj  BPTREE_OPLIST2(dbg_adj, M_BASIC_OPLIST, M_OPL_array_core_locs)

//BPTREE_DEF2(debruin, 4, hashrb_t, M_OPEXTEND(M_POD_OPLIST, HASH(hashrb_rehash2), EQUAL(hashrb_eq), CMP(hashrb_cmp)), dbg_adj_t, M_OPL_dbg_adj)
BPTREE_DEF2(debruin, 4, struct lcmer_t, M_POD_OPLIST, dbg_adj_t, M_OPL_dbg_adj)



#define SHIFT_COUNT (64/BUFFER_SIZE)
uint64_t hash_move_right(uint64_t h, uint64_t n){
	return h ^ n;	
}
uint64_t hash_move_left(uint64_t h, uint64_t n){
	return h ^ n;
}


// Maybe a bloom filter 
#define MIN_REQUIRED_CORE_COUNT 2
BPTREE_DEF2(counter, 4, uint64_t, size_t)

BPTREE_DEF2(lcmer_counter, 4, struct lcmer_t, M_POD_OPLIST, size_t, M_BASIC_OPLIST)


void lcmer_shift_left(struct lcmer_t *l){
	for(int i = 0; i < BUFFER_SIZE -1; ++i){
		l->data[i]=l->data[i+1];
	}
}
void dfs_step(debruin_t *r, dbg_adj_t *adj, lcmer_counter_t *c, struct lcmer_t *lc, FILE *out, uint64_t prev_hash){

	size_t *visited = lcmer_counter_safe_get(*c, *lc);
	if(*visited == 0){
		++(*visited);
		dbg_adj_it_t ait;
		lcmer_shift_left(lc);
		for(dbg_adj_it(ait, *adj); !dbg_adj_end_p(ait); dbg_adj_next(ait)){
			uint64_t new_core = *dbg_adj_ref(ait)->key_ptr;
			lc->data[BUFFER_SIZE-1] = new_core;
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

void lcp_space_graph_construct_from_genome(struct ref_seq *seqs, FILE *out) {

//	hashrb_t rb = {0};
	
	counter_t counter_dict;
	counter_init(counter_dict);

	debruin_t left_dbg;
	debruin_init(left_dbg);

	debruin_t right_dbg;
	debruin_init(right_dbg);

	for(int ch_i = 0; ch_i < seqs->size; ++ch_i){
		int core_i;
		for(core_i = 0; core_i < seqs->chrs[ch_i].cores_size; ++core_i){
			size_t *count = counter_safe_get( counter_dict, seqs->chrs[ch_i].cores[core_i].id);
			++(*count);
		}

	}

	for(int ch_i = 0; ch_i < seqs->size; ++ch_i){
		int core_i = seqs->chrs[ch_i].cores_size;
		int fill = 0;
		struct lcmer_t lcl = {0};
		
		//while(core_i < seqs->chrs[ch_i].cores_size && fill < BUFFER_SIZE){
		if(0) { // reversed
			while(core_i >= 0 && fill < BUFFER_SIZE){
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
			while(core_i < seqs->chrs[ch_i].cores_size && fill < BUFFER_SIZE){
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
		for(int i = 0; i < BUFFER_SIZE/2;++i){
			int tmp = lcl.data[i];
			lcl.data[i] = lcl.data[BUFFER_SIZE-1-i];
			lcl.data[BUFFER_SIZE-1-i] = tmp;
		}
		uint64_t prev = 0;
		if(0){//reversed
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
				array_core_locs_push_back(*dbg_adj_safe_get(*rad, core_in), cl);

				if(dbg_adj_size(*rad) == 1 && array_core_locs_size(*dbg_adj_safe_get(*rad, core_in))== 1){
					fprintf(out, "S\t");
					fprintf(out, "%lu", hv);
					for(int i = 0; i < BUFFER_SIZE; ++i){
						fprintf(out, "-%lu", lcl.data[i]);
					}
					fprintf(out, "\t*\tRC:i:%d\n", ch_i);
				}
				if(prev!=0){
					fprintf(out, "L\t%lu\t+\t%lu\t+\t*\tID:Z:%d\n", prev, hv, ch_i);
				}

				for(int i = 0; i < BUFFER_SIZE -1; ++i){
					lcl.data[i]=lcl.data[i+1];
				}
				lcl.data[BUFFER_SIZE-1] = core_in;
				prev = hv;
			}
		}
	}

	lcmer_counter_t cd;
	lcmer_counter_init(cd);
	dfs(&right_dbg, &cd, out);

	counter_clear(counter_dict);
	debruin_clear(right_dbg);
	debruin_clear(left_dbg);
}

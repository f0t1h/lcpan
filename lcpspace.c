#include "lcpspace.h"
#include "ringbuffer.h"
#include "fa_parser.h"
#include "m-bptree.h"
#include "m-array.h"

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
    case 2: k ^= (tail[1] << 8);

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
#define BUFFER_SIZE 8
DEFINE_HASHED_RING_BUFFER(hashrb, HASH_COUNT, BUFFER_SIZE, FNVHashInt, MurmurHash3_int) 

struct core_loc {
	int start;
	int end;
	uint64_t tid;
};


ARRAY_DEF(array_core_locs, struct core_loc, M_POD_OPLIST)
#define M_OPL_array_core_locs  ARRAY_OPLIST(array_core_locs, M_POD_OPLIST)
BPTREE_DEF2(dbg_adj, 4, uint64_t, M_BASIC_OPLIST, array_core_locs_t, M_OPL_array_core_locs)
#define M_OPL_dbg_adj  BPTREE_OPLIST2(dbg_adj, M_BASIC_OPLIST, M_OPL_array_core_locs)
BPTREE_DEF2(debruin, 4, uint64_t, M_BASIC_OPLIST,dbg_adj_t, M_OPL_dbg_adj)


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

void dfs(debruin_t *l, debruin_t *r){

	debruin_it_t it;
	for(debruin_it(it, *r); !debruin_end_p(it); debruin_next(it)){
		
	}
}

void lcp_space_graph_construct_from_genome(struct ref_seq *seqs, FILE *out) {

	hashrb rb = {0};
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
		int core_i = 0;
		int fill = 0;
		while(fill < BUFFER_SIZE){
			size_t *count = counter_safe_get( counter_dict, seqs->chrs[ch_i].cores[core_i].id);
			if(*count >= MIN_REQUIRED_CORE_COUNT){
				++fill;
				hashrb_put( &rb, seqs->chrs[ch_i].cores[core_i].id);
			}
			++core_i;
		}

		for(; core_i < seqs->chrs[ch_i].cores_size; ++core_i){
			size_t *count = counter_safe_get( counter_dict, seqs->chrs[ch_i].cores[core_i].id);
			if(*count >= MIN_REQUIRED_CORE_COUNT){

				uint64_t core_out = 0;
				hashrb_get( &rb, &core_out);
				uint64_t hv = hashrb_rehash( &rb); 
				uint64_t core_in = seqs->chrs[ch_i].cores[core_i].id;
				dbg_adj_t *lad = debruin_safe_get(left_dbg, hv);
				struct simple_core sc = seqs->chrs[ch_i].cores[core_i];
				struct core_loc cl = {sc.start, sc.end, ch_i};
				array_core_locs_push_back(*dbg_adj_safe_get(*lad, core_out), cl);

				dbg_adj_t *rad = debruin_safe_get(right_dbg, hv);
				array_core_locs_push_back(*dbg_adj_safe_get(*rad, core_in), cl);

				hashrb_put( &rb, core_in);
			}
		}
	}


	counter_clear(counter_dict);
	debruin_clear(right_dbg);
	debruin_clear(left_dbg);
}

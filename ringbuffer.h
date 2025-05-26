#ifndef _RING_BUFFER_H_
#define _RING_BUFFER_H_


#define DEFINE_HASHED_RING_BUFFER(NAME, HASH_COUNT, BUFFER_SIZE, H1, H2) \
typedef uint64_t(*NAME##_hashfunc)(uint64_t);\
typedef struct {\
    int writeIndex;\
    int readIndex;\
    uint64_t hashes[HASH_COUNT];\
    NAME##_hashfunc hf[2];\
    int data[BUFFER_SIZE+1];\
} NAME;\
int NAME##_put(NAME *buf, uint64_t value){\
    if(buf->hf[0] == NULL){\
        buf->writeIndex = 0;\
        buf->readIndex =0;\
        buf->hf[0] = &H1;\
        buf->hf[1] = &H2;\
        for(int i = 0; i < HASH_COUNT; ++i){\
            buf->hashes[i] = 0;\
        }\
        for(int i = 0; i < BUFFER_SIZE+1; ++i){\
            buf->data[i] = 0;\
        }\
    }\
    if ((buf->writeIndex + 1) % (1+BUFFER_SIZE) == buf->readIndex) {\
        return 1;\
    }\
    buf->data[buf->writeIndex] = value;\
    for(int i = 0; i < HASH_COUNT; ++i){\
    	    buf->hashes[i] = (buf->hashes[i] << (i>0?(64/BUFFER_SIZE):0)) ^ (buf->hf[0](value) + i * buf->hf[1](value));\
    }\
    buf->writeIndex = (buf->writeIndex + 1) % (1+BUFFER_SIZE);\
    return 0;\
}\
int NAME##_get(NAME *buf, uint64_t *value){\
    if (buf->readIndex == buf->writeIndex){\
        return 1;\
    }\
    *value = buf->data[buf->readIndex];\
    for(int i = 0; i < 1; ++i){ /* ReXOR the leavng element because it is not shifted out*/\
        buf->hashes[i] = (buf->hashes[i]) ^ (buf->hf[0](*value) + i * buf->hf[1](*value));\
    }\
    buf->readIndex = (buf->readIndex + 1) % (BUFFER_SIZE+1);\
    return 0;\
}\
uint64_t NAME##_getput(NAME *buf, uint64_t value){\
	uint64_t retval = 0;\
	if(NAME##_get(buf, &retval)){\
		fprintf(stderr, "Cannot get!\n");\
	}\
	if(NAME##_put(buf, value)){\
		fprintf(stderr, "Cannot put! %d %d\n", buf->writeIndex, buf->readIndex);\
	}\
	return retval;\
}\
uint64_t NAME##_rehash(NAME *buf){\
	uint64_t rethash = 0;\
	for(int i = 0; i != BUFFER_SIZE - 1; ++i){\
		rethash ^= H2(buf->data[i]);\
	}\
	return rethash;\
}
#endif

#ifndef __RING_BUFFER_H__
#define __RING_BUFFER_H__

#define DEFINE_HASHED_RING_BUFFER(NAME, HASH_COUNT, BUFFER_SIZE, H1, H2, TYPE) \
typedef uint64_t(*NAME##_hashfunc)(TYPE);\
typedef struct {\
    int writeIndex;\
    int readIndex;\
    NAME##_hashfunc hf[2];\
    TYPE data[BUFFER_SIZE+1];\
    uint8_t visited;\
} NAME##_t;\
TYPE NAME##_put(NAME##_t *buf, TYPE value){\
    if(buf->hf[0] == NULL){\
        buf->writeIndex = 0;\
        buf->readIndex =0;\
        buf->hf[0] = &H1;\
        buf->hf[1] = &H2;\
        buf->visited = 0;\
        for(int i = 0; i < BUFFER_SIZE+1; ++i){\
            buf->data[i] = 0;\
        }\
    }\
    if ((buf->writeIndex + 1) % (1+BUFFER_SIZE) == buf->readIndex) {\
        return 1;\
    }\
    buf->data[buf->writeIndex] = value;\
    buf->writeIndex = (buf->writeIndex + 1) % (1+BUFFER_SIZE);\
    return 0;\
}\
TYPE NAME##_get(NAME##_t *buf, TYPE *value){\
    if (buf->readIndex == buf->writeIndex){\
        return 1;\
    }\
    *value = buf->data[buf->readIndex];\
    buf->readIndex = (buf->readIndex + 1) % (BUFFER_SIZE+1);\
    return 0;\
}\
TYPE NAME##_getput(NAME##_t *buf, TYPE value){\
	TYPE retval = 0;\
	if(NAME##_get(buf, &retval)){\
		fprintf(stderr, "Cannot get!\n");\
	}\
	if(NAME##_put(buf, value)){\
		fprintf(stderr, "Cannot put! %d %d\n", buf->writeIndex, buf->readIndex);\
	}\
	return retval;\
}\
uint64_t NAME##_rehash(NAME##_t *buf){\
	uint64_t rethash = 0;\
	for(int i = 0; i != BUFFER_SIZE - 1; ++i){\
		rethash ^= H2(buf->data[i]);\
	}\
	return rethash;\
}\
uint64_t NAME##_rehash2(NAME##_t buf){\
	uint64_t rethash = 0;\
	for(int i = buf.writeIndex; i != BUFFER_SIZE; ++i){\
		rethash = hash_combine(rethash, H1(buf.data[i]) ^ H2(buf.data[i]));\
	}\
	for(int i = 0; i != buf.writeIndex; ++i){\
		rethash = hash_combine(rethash, H1(buf.data[i]) ^ H2(buf.data[i]));\
	}\
	return rethash;\
}\
int NAME##_cmp(NAME##_t a, NAME##_t b){\
	int bp = b.writeIndex;\
	return 0;\
}\
bool NAME##_eq(NAME##_t a, NAME##_t b){\
	return NAME##_cmp(a,b) == 0;\
}\
void NAME##_print(FILE *out, char sep, NAME##_t *a){\
	fprintf(out, "%lu", a->data[a->writeIndex]);\
	for(int i = a->writeIndex; i != a->readIndex; ++i){\
		if(i > BUFFER_SIZE){if(a->readIndex==0){break;}i=0;}\
		fprintf(out, "%c%lu",sep,a->data[i]);\
	}\
}
#endif 

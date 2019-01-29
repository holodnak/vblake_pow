#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef unsigned __int64 uint64_t;

// The modified for bias-free initialization vector for the blake2b algorithm
uint64_t vBlake_iv[] =  {
	0x4BBF42C1F006AD9DULL, 0x5D11A8C3B5AEB12EULL,
	0xA64AB78DC2774652ULL, 0xC67595724658F253ULL,
	0xB8864E79CB891E56ULL, 0x12ED593E29FB41A1ULL,
	0xB1DA3AB63C60BAA8ULL, 0x6D20E50C1F954DEDULL
};

// Sigma array for use in G mixing function, note re-addition of
// original 4 last sigma lines from BLAKE to allow 16 total rounds
uint8_t sigma[16][16] = {
	{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
	{ 14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3 },
	{ 11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4 },
	{ 7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8 },
	{ 9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13 },
	{ 2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9 },
	{ 12, 5, 1, 15, 14, 13, 4, 10, 0, 7, 6, 3, 9, 2, 8, 11 },
	{ 13, 11, 7, 14, 12, 1, 3, 9, 5, 0, 15, 4, 8, 6, 2, 10 },
	{ 6, 15, 14, 9, 11, 3, 0, 8, 12, 2, 13, 7, 1, 4, 10, 5 },
	{ 10, 2, 8, 4, 7, 6, 1, 5, 15, 11, 9, 14, 3, 12, 13, 0 },
	{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
	{ 14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3 },
	{ 11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4 },
	{ 7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8 },
	{ 9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13 },
	{ 2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9 }
};

// The re-introduced constants, modified to also be bias-free
uint64_t vBlake_c[] = {
	0xA51B6A89D489E800ULL, 0xD35B2E0E0B723800ULL,
	0xA47B39A2AE9F9000ULL, 0x0C0EFA33E77E6488ULL,
	0x4F452FEC309911EBULL, 0x3CFCC66F74E1022CULL,
	0x4606AD364DC879DDULL, 0xBBA055B53D47C800ULL,
	0x531655D90C59EB1BULL, 0xD1A00BA6DAE5B800ULL,
	0x2FE452DA9632463EULL, 0x98A7B5496226F800ULL,
	0xBAFCD004F92CA000ULL, 0x64A39957839525E7ULL,
	0xD859E6F081AAE000ULL, 0x63D980597B560E6BULL
};

void compress(uint64_t *h, uint8_t *b);
void recombineB2Bh(uint8_t *output, uint64_t *h);
uint64_t bytesToLong(uint8_t *bytes);
void B2B_G(uint64_t *v, int a, int b, int c, int d, uint64_t x, uint64_t y, uint64_t c1, uint64_t c2);

int hash(uint8_t *output, uint8_t *input, int inputlen) {
	// BEGIN INIT
	uint64_t h[8];
	uint8_t b[64];

	h[0] = vBlake_iv[0];
	h[1] = vBlake_iv[1];
	h[2] = vBlake_iv[2];
	h[3] = vBlake_iv[3];
	h[4] = vBlake_iv[4];
	h[5] = vBlake_iv[5];
	h[6] = vBlake_iv[6];
	h[7] = vBlake_iv[7];

	// outlen = 24, as VeriBlock uses a 192-bit hash
	h[0] ^= (uint64_t)(0x01010000 ^ 0x18);

	// END INIT
	// BEGIN UPDATE

	for (int i = 0; i < inputlen; i++) {
		b[i] = input[i];
	}

	compress(h, b);
	// END UPDATE

	recombineB2Bh(output, h);
	return 0;
}

void compress(uint64_t *h, uint8_t *b) {
	uint64_t v[16];
	uint64_t m[16];

	for (int i = 0; i < 8; i++) {
		v[i] = h[i];
		v[i + 8] = vBlake_iv[i];
	}

	v[12] ^= 64; // Input count low
	v[13] ^= 0;  // Input count high (no overflow, therefore 0)

	v[14] ^= (long)(-1); // f[0] = 0xFF..FF
	v[15] ^= 0;          // f[1] = 0x00..00

	memset(m, 0, sizeof(uint64_t) * 16);
	for (int i = 0; i < 8; i++) {
		m[i] = bytesToLong(&b[i * 8]);
	}

	// Using 16 rounds of the Blake2 G function, drawing on the additional 4 rows
	// of sigma from reference BLAKE implementation
	for (int i = 0; i < 16; i++) {
		B2B_G(v, 0, 4, 8, 12, m[sigma[i][1]], m[sigma[i][0]], vBlake_c[sigma[i][1]], vBlake_c[sigma[i][0]]);
		B2B_G(v, 1, 5, 9, 13, m[sigma[i][3]], m[sigma[i][2]], vBlake_c[sigma[i][3]], vBlake_c[sigma[i][2]]);
		B2B_G(v, 2, 6, 10, 14, m[sigma[i][5]], m[sigma[i][4]], vBlake_c[sigma[i][5]], vBlake_c[sigma[i][4]]);
		B2B_G(v, 3, 7, 11, 15, m[sigma[i][7]], m[sigma[i][6]], vBlake_c[sigma[i][7]], vBlake_c[sigma[i][6]]);
		B2B_G(v, 0, 5, 10, 15, m[sigma[i][9]], m[sigma[i][8]], vBlake_c[sigma[i][9]], vBlake_c[sigma[i][8]]);
		B2B_G(v, 1, 6, 11, 12, m[sigma[i][11]], m[sigma[i][10]], vBlake_c[sigma[i][11]], vBlake_c[sigma[i][10]]);
		B2B_G(v, 2, 7, 8, 13, m[sigma[i][13]], m[sigma[i][12]], vBlake_c[sigma[i][13]], vBlake_c[sigma[i][12]]);
		B2B_G(v, 3, 4, 9, 14, m[sigma[i][15]], m[sigma[i][14]], vBlake_c[sigma[i][15]], vBlake_c[sigma[i][14]]);
	}

	// Update h[0 .. 7]
	for (int i = 0; i < 8; i++) {
		h[i] ^= v[i] ^ v[i + 8];
	}

	h[0] ^= h[3] ^ h[6];
	h[1] ^= h[4] ^ h[7];
	h[2] ^= h[5];
}

/**
* Rotate x, a 64-bit long, to the right by y places
*/
uint64_t ROTR64(uint64_t x, int y) {
	return (x >> y) | (x << (64 - y));
}

/**
* The G Mixing function from the Blake2 specification.
*/
void B2B_G(uint64_t *v, int a, int b, int c, int d, uint64_t x, uint64_t y, uint64_t c1, uint64_t c2) {

	v[a] = v[a] + v[b] + (x ^ c1);
	v[d] ^= v[a];
	v[d] = ROTR64(v[d], 60);
	v[c] = v[c] + v[d];
	v[b] = ROTR64(v[b] ^ v[c], 43);
	v[a] = v[a] + v[b] + (y ^ c2);
	v[d] = ROTR64(v[d] ^ v[a], 5);
	v[c] = v[c] + v[d];
	v[b] = ROTR64(v[b] ^ v[c], 18);

	// X'Y'Z' + X'YZ + XY'Z + XYZ'    LUT: 10010110
	v[d] ^= (~v[a] & ~v[b] & ~v[c]) | (~v[a] & v[b] & v[c]) |
		(v[a] & ~v[b] & v[c]) | (v[a] & v[b] & ~v[c]);

	// X'Y'Z + X'YZ' + XY'Z' + XYZ    LUT: 01101001
	v[d] ^= (~v[a] & ~v[b] & v[c]) | (~v[a] & v[b] & ~v[c]) |
		(v[a] & ~v[b] & ~v[c]) | (v[a] & v[b] & v[c]);
}

/**
* Convert 8 bytes in a byte[] to a single long, maintaining endianness.
*/
uint64_t bytesToLong(uint8_t *bytes) {
	uint64_t result =
		((((uint64_t)bytes[7]) & 0xFF) << 56) |
		((((uint64_t)bytes[6]) & 0xFF) << 48) |
		((((uint64_t)bytes[5]) & 0xFF) << 40) |
		((((uint64_t)bytes[4]) & 0xFF) << 32) |
		((((uint64_t)bytes[3]) & 0xFF) << 24) |
		((((uint64_t)bytes[2]) & 0xFF) << 16) |
		((((uint64_t)bytes[1]) & 0xFF) << 8) |
		((((uint64_t)bytes[0]) & 0xFF));
	return result;
}

/**
* Changes the h[0 .. 2] array back to a byte[0 .. 23] array
*  for the final hash output.
*/
void recombineB2Bh(uint8_t *output, uint64_t *h) {
	for (int i = 0; i < 3; i++) {
		output[i * 8 + 0] = (uint8_t)(h[i] >> 0);
		output[i * 8 + 1] = (uint8_t)(h[i] >> 8);
		output[i * 8 + 2] = (uint8_t)(h[i] >> 16);
		output[i * 8 + 3] = (uint8_t)(h[i] >> 24);
		output[i * 8 + 4] = (uint8_t)(h[i] >> 32);
		output[i * 8 + 5] = (uint8_t)(h[i] >> 40);
		output[i * 8 + 6] = (uint8_t)(h[i] >> 48);
		output[i * 8 + 7] = (uint8_t)(h[i] >> 56);
	}
}

void printData(unsigned char *data, int size)
{
	int i;
	for (i = 0; i<size; i++)
	{
		printf("%02X", data[i]);
		if ((i + 1) % 16 == 0) printf("\n");
		else if ((i + 1) % 8 == 0) printf(" - ");
		else if ((i + 1) % 4 == 0) printf(" ");
	}
	printf("\n");
}

void bswap(unsigned char *b, int len)
{
	while (len) {
		unsigned char t[4], i;

		for (i = 0; i < 4; i++)	t[i] = b[i];
		for (i = 0; i < 4; i++)	b[i] = t[3 - i];
		b += 4;
		len -= 4;
	}
}

int main(int ac, char **av)
{
	uint32_t data_in[16] = {
		0x0005A619, 0x0001FFC3, 0xAFFCDBD5, 0xC2764320,
		0xB38A0D92, 0x562B8AE6, 0xFCAEA270, 0xE66D7E22,
		0xDC07A7D3, 0xF00C03A5, 0xB6CE5F6B, 0x1CDD8F8D,
		0x50D72757, 0x5C4ACF75, 0x063A3026, 0xc471f45f
	};
	uint32_t expected_result[2] = { 0x00000000, 0x76A0569C };
	//000000009C56A076

	uint64_t h[8];
	uint8_t in[64];

	memcpy(in, data_in, 64);
	bswap(in, 64);

	hash((uint8_t*)h, in, 64);

	printf("input 64-byte data:\n");
	printData((uint8_t*)in, 64);

	printf("final 24-byte hash:\n");
	printData((uint8_t*)h, 24);

	if (memcmp(h, expected_result, 8) == 0) {
		printf("\nhash success!\n");
	}

	else {
		printf("hash generation failed.\n");
	}

	system("pause");
	return 0;
}
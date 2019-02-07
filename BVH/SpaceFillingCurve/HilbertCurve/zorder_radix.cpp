#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <bitset>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>



static uint32_t expandBits(unsigned int v)
{
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

uint32_t morton3d(float x, float y, float z)
{
    float norm =  sqrt(x*x + y*y + z*z);
    x /= norm; y /= norm; z /= norm;

    x = std::min(std::max(x * 1024.0f, 0.0f), 1023.0f);
    y = std::min(std::max(y * 1024.0f, 0.0f), 1023.0f);
    z = std::min(std::max(z * 1024.0f, 0.0f), 1023.0f);
    uint32_t xx = expandBits((uint32_t)x);
    uint32_t yy = expandBits((uint32_t)y);
    uint32_t zz = expandBits((uint32_t)z);
    return xx * 4 + yy * 2 + zz;
}

// replaced byte with bitsOffset to avoid *8 operation in loop
static void radix(short byteOffset, const unsigned int N, uint32_t* source, uint32_t* dest)
{
	// suppressed the need for index as it is reported in count
	uint32_t count[256];
	// added temp variables to simplify writing, understanding and compiler
	// optimization job most of them will be allocated as registers
	uint32_t *sp, *cp, s, c, i;
	uint8_t *bp;

	// faster than MemSet
	cp = count;
	for (i = 256; i > 0; --i, ++cp)
		*cp = 0;

	// count occurences of every byte value
	bp = ((uint8_t*) source) + byteOffset;
	for (i = N; i > 0; --i, bp += 4) {
		cp = count + *bp;
		++(*cp);
	}

	// transform count into index by summing elements and storing into same array
	s = 0;
	cp = count;
	for (i = 256; i > 0; --i, ++cp) {
		c = *cp;
		*cp = s;
		s += c;
	}

	// fill dest with the right values in the right place
	bp = ((uint8_t*) source) + byteOffset;
	sp = source;
	for (i = N; i > 0; --i, bp += 4, ++sp) {
		cp = count + *bp;
		dest[*cp] = *sp;
		++(*cp);
	}
}

//static void radix_sort (uint32_t *source, unsigned N)
void radix_sort(uint32_t* source, const unsigned int N)
{
	// allocate heap memory to avoid the need of additional parameter
	uint32_t* temp = (uint32_t*) malloc(N*sizeof(uint32_t));
	assert (temp != NULL);

	radix (0, N, source, temp);
	radix (1, N, temp, source);
	radix (2, N, source, temp);
	radix (3, N, temp, source);

	free (temp);
}

//static void check_order (uint32_t *data, unsigned N)
void check_order(uint32_t* data, const unsigned int N)
{
    unsigned int n = N;
	// only signal errors if any (should not be)
	for( int i = 1; i < n; i++ )
        assert( data[i-1] <= data[i] );
    /*for (--n ; n > 0; --n, ++data)
		assert (data[0] <= data[1]);*/
}

//return index of duplicate
int check_duplicates(uint32_t* data, const unsigned int N)
{
    unsigned int n = N;
    for( int i = 1; i < n; i++ )
    {
        if( data[i-1] == data[i] )
            return i;
    }
    return 0;
}//can combine these into since uniqueness check function


static void mkdirTree(std::string sub, std::string dir)
{
    if(sub.length() == 0)
        return;

    int i = 0;
    for( i; i < sub.length(); i++)
    {
        dir += sub[i];
        if (sub[i] == '/')
            break;
    }

    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if(i+1 < sub.length())
        mkdirTree(sub.substr(i+1), dir);
}


void create_directory(std::string new_dir)
{
    struct stat st;
    int status = stat(new_dir.c_str(), &st);
    if( status != 0 && !S_ISDIR(st.st_mode) )
        mkdirTree(new_dir, "");
    /*else
        std::cout << "Directory already exists.\n";*/
}





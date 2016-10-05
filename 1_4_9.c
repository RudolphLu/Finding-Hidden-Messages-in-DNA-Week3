/*
 * Given a collection of strings Dna and an integer d, a k-mer is a (k,d)-motif 
 * if it appears in every string from Dna with at most d mismatches. For example, 
 * the implanted 15-mer in the strings above represents a (15,4)-motif.
 *
 * Delete the latest "\n" from input file
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEBUG

typedef unsigned char    UINT8;
typedef unsigned int     UINT32;


#define MAX_ARRAY_SIZE   1000000

#define SCAN_KMER        0
#define SCAN_DNA         1


UINT8  DNA_SRC[MAX_ARRAY_SIZE];
UINT32 NUMBER_OF_DNAS;
UINT32 SIZE_OF_DNA = 0;
UINT32 DNA_TOTAL_LENGTH;
UINT32 KMER_LENGTH = 0;
UINT32 HAMMING_DISTANCE = 0;

UINT32 HASH_TABLE_LENGTH = 1;
UINT32 *pALL_KEY;
UINT32 alpha[256] = {0};

typedef struct pattern{
    UINT32 kmer;
    UINT32 key;
}pattern_t;

typedef struct dna{
    UINT32 number;
    UINT32 length;
    UINT8  *src;
}dna_t;


/*
 * Key base.
 */
UINT32 hamming_distance_bykey(UINT32 key1,UINT32 key2,UINT32 len)
{
    UINT32 distance = 0, i = 0;
    for(i=0;i<len;i++){
        if((key1 & (0x3<<(i*2)))!=(key2 & (0x3<<(i*2))))
            distance++; 
    }
    return distance;
}


UINT32 gen_hash_key(UINT8 *p,UINT32 len)
{
    UINT32 key = 0;
    UINT32 i=0;
    UINT32 v;
    
    for(i=0;i<len;i++){        
        v = (UINT32)alpha[p[i]];        
        key |= (v) << (2*i);
    }
    return key;
}

/*
 * Return the min sum of hamming distance in all DNA strings.
 */
UINT32 distance_between_pattern_and_strings(pattern_t *p,dna_t *d)
{
    UINT32 i,j;
    UINT32 sum_hamming_distance=0;
    UINT32 key_dna,key_pattern = p->key;
    UINT32 distance = 0xffffffff,key_distance;
    UINT32 kmer = p->kmer;
    UINT8  *src = d->src;
    
    
    for(i=0,distance=0xffffffff;i<d->number;i++,distance=0xffffffff){
        for(j=0;j<=(SIZE_OF_DNA-kmer);j++){
            //compare key
            key_dna = gen_hash_key(src+(i*d->length+j),kmer);
            key_distance = hamming_distance_bykey(key_pattern,key_dna,kmer);
            if(distance>key_distance)
                distance=key_distance;
        }
        sum_hamming_distance += distance;
    }
    
    return sum_hamming_distance;
}


void read_test_data (const char* file_name,UINT8 *src ,int maxsize) 
{    
    UINT32 phase = SCAN_KMER,i,src_cnt=0;
    UINT8 c;    
    FILE* file = fopen (file_name, "r");
    
    if (!file) return;

    HAMMING_DISTANCE = KMER_LENGTH = 0;
    HASH_TABLE_LENGTH = 1;
    
    do{
        c = fgetc(file);    
        if( feof(file) ){
           break ;
        }
        
        if((0xa==c)||(0xd==c)){
            phase = SCAN_DNA;
            if(0xa==c){
              NUMBER_OF_DNAS++;
              SIZE_OF_DNA = 0;
            }
            continue;
        }

        if(SCAN_DNA == phase){
            if( (c == 'G') || (c == 'C') || (c=='A') || (c=='T')){
                src[src_cnt++] = c;
                SIZE_OF_DNA ++; // Each size of DNAS should be equal.
            }
        }else if(SCAN_KMER == phase){ 
            KMER_LENGTH *= 10;          
            KMER_LENGTH += c-'0';
        }
                
    }while(1);

    fclose (file);

    for(i=0;i<KMER_LENGTH;i++)
        HASH_TABLE_LENGTH *= 4;    

    DNA_TOTAL_LENGTH = src_cnt;
}

void printpattern(UINT32 key)
{
    UINT32 l;
    for(l=0;l<KMER_LENGTH;l++){
        if((key&0x3) == 0x0) printf("A");
        else if((key&0x3) == 0x1) printf("C");
        else if((key&0x3) == 0x2) printf("G");
        else if((key&0x3) == 0x3) printf("T");
        key>>=2;
    }
    printf(" ");
}

void main()
{
    UINT32 i,min_distance,sum_of_min_d;
    UINT32 match_cnt = 0;
    dna_t  dna_src;
    pattern_t cmp_pattern;
    read_test_data("dataset_158_9.txt",DNA_SRC,MAX_ARRAY_SIZE);
    
    dna_src.number = NUMBER_OF_DNAS;
    dna_src.src = DNA_SRC;
    dna_src.length = SIZE_OF_DNA;
    
#ifdef DEBUG
    printf("num:%d sz:%d kmer_sz:%d hd:%d \n",NUMBER_OF_DNAS,SIZE_OF_DNA,KMER_LENGTH,HAMMING_DISTANCE);
#endif
    // Init hash index.
    alpha['A'] = 0;
    alpha['C'] = 1;
    alpha['G'] = 2;
    alpha['T'] = 3;
    
    //Init hash table.
    pALL_KEY = malloc(HASH_TABLE_LENGTH);
    memset(pALL_KEY,0,HASH_TABLE_LENGTH*sizeof(UINT32));
    
    for(i=0,min_distance=0xffffffff;i<HASH_TABLE_LENGTH;i++){
        cmp_pattern.key = i;
        cmp_pattern.kmer = KMER_LENGTH;
        sum_of_min_d = distance_between_pattern_and_strings(&cmp_pattern,&dna_src);
        pALL_KEY[cmp_pattern.key] = sum_of_min_d;
        if(min_distance>sum_of_min_d)
            min_distance = sum_of_min_d;
    }


    for(i=0;i<HASH_TABLE_LENGTH;i++){
        if(min_distance==pALL_KEY[i])
            printpattern(i);
    }
}
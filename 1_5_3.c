/*
 * Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string
 *
 * Delete the latest "\n" from input file
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define DEBUG

typedef unsigned char    UINT8;
typedef unsigned int     UINT32;


#define MAX_ARRAY_SIZE   1000000

#define SCAN_DNA         0
#define SCAN_KMER        1
#define SCAN_PROFILE     2



UINT8  DNA_SRC[MAX_ARRAY_SIZE];
double PROFILE_PATTERN[MAX_ARRAY_SIZE];
UINT32 NUMBER_OF_PATTERNS;
UINT32 SIZE_OF_PATTERN = 0;
UINT32 DNA_TOTAL_LENGTH;
UINT32 KMER_LENGTH = 0;
UINT32 alpha[256] = {0};


void read_test_data (const char* file_name,UINT8 *src ,double *profile,int maxsize) 
{    
    UINT32 phase = SCAN_DNA,count=0;
    double i;
    UINT8 c; 
    UINT32 j;   
    FILE* file = fopen (file_name, "r");
    
    if (!file) return;

    SIZE_OF_PATTERN = KMER_LENGTH = DNA_TOTAL_LENGTH = 0;
    
    do {
        if(SCAN_DNA == phase){
            c = fgetc(file);
            if( (c == 'G') || (c == 'C') || (c=='A') || (c=='T')){
                src[DNA_TOTAL_LENGTH ++] = c;                
            }else if(c == 0xd){
                continue;
            }else if(c == 0xa){
                phase = SCAN_KMER;
            }
        }
        
        if(SCAN_KMER == phase){
            c = fgetc(file);
            if(c == 0xd){
                continue;
            }else if(c == 0xa){
                phase = SCAN_PROFILE;
                continue;
            }
            KMER_LENGTH *= 10;
            KMER_LENGTH += c-'0';
        }
        
        if(SCAN_PROFILE == phase){
           fscanf (file, "%lf", &i);
           profile[SIZE_OF_PATTERN++] = i;
        }
        
    } while (!feof (file) && count < maxsize);

    NUMBER_OF_PATTERNS = SIZE_OF_PATTERN/KMER_LENGTH;
    fclose (file);
}

UINT32 caculate_pr(UINT8 *dna)
{
    UINT32 i;
    double multiple = 1;
    
    for(i=0,multiple=1;i<KMER_LENGTH;i++){             
        multiple*=PROFILE_PATTERN[alpha[dna[i]]*KMER_LENGTH+i];
    }
    
    return multiple;
}


void main()
{
    UINT32 i=0;
    double pr_multi,max_pr=0;
    UINT32 max_pr_position;
    read_test_data("dataset_158_9.txt",DNA_SRC,PROFILE_PATTERN,MAX_ARRAY_SIZE);
#ifdef DEBUG
    printf("num_p:%d sz_p:%d dna_sz:%d kmer:%d \n",NUMBER_OF_PATTERNS,SIZE_OF_PATTERN,DNA_TOTAL_LENGTH,KMER_LENGTH);
    for(i=0;i<SIZE_OF_PATTERN;i++)
        printf("%d ",PROFILE_PATTERN[i]);
#endif
    // Init hash index.
    alpha['A'] = 0;
    alpha['C'] = 1;
    alpha['G'] = 2;
    alpha['T'] = 3;
    
    
    for(i=0;i<=(DNA_TOTAL_LENGTH-KMER_LENGTH);i++){
        pr_multi=caculate_pr(DNA_SRC+i);
        if(max_pr<pr_multi){
            max_pr = pr_multi;
            max_pr_position = i;
        }            
    }
    
    for(i=max_pr_position;i<(max_pr_position+KMER_LENGTH);i++){
        printf("%c",DNA_SRC[i]);
    }
    
}
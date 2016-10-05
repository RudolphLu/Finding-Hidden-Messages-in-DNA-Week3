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
#define SCAN_HAMMING_D   1
#define SCAN_DNA         2

#define NEW_PATTERN      0xffffffff

UINT8  DNA_SRC[MAX_ARRAY_SIZE];
UINT32 NUMBER_OF_DNAS;
UINT32 SIZE_OF_DNA = 0;
UINT32 DNA_TOTAL_LENGTH;
UINT32 KMER_LENGTH = 0;
UINT32 HAMMING_DISTANCE = 0;

UINT32 HASH_TABLE_LENGTH = 1;
UINT32 *pKD_MOTIF;
UINT32 alpha[256] = {0};

UINT32 gen_hash_key(unsigned char *p,int len)
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

        if(' ' == c){
            phase = SCAN_HAMMING_D;
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
        }else if(SCAN_HAMMING_D == phase){
            HAMMING_DISTANCE *= 10;
            HAMMING_DISTANCE += c-'0';
        }
                
    }while(1);

    fclose (file);

    for(i=0;i<KMER_LENGTH;i++)
        HASH_TABLE_LENGTH *= 4;    

    DNA_TOTAL_LENGTH = src_cnt;
}

/*
 * Generate (k,d)-motif  
 * offset: Current position of nucleotides in pattern.
 * hd:     Hamming distance.
 * len:    The k-mer.
 */
void gen_kd_motif(UINT32 key,UINT32 offset,UINT32 hd,UINT32 len)
{
    UINT32 i = 0;
    UINT32 nucleotide = (key>>(2*offset)) & 0x3;
    UINT32 mask = 0x3<<(2*offset);
    UINT32 temp_key = key;
    
    
    if(offset==len)
        return ;

    if(hd!=0){
        for(i=0;i<=3;i++){
            if(i!=nucleotide){
                 temp_key &= ~mask;
                 temp_key |= i<<(2*offset);
                 
                 //Add new key.
                 if(0 == pKD_MOTIF[temp_key])
                     pKD_MOTIF[temp_key] = NEW_PATTERN;

                 gen_kd_motif(temp_key,offset+1,hd-1,len);
            }
        }                    
        gen_kd_motif(key,offset+1,hd,len);
    }else
        return;

}

/*
 * Find out the (k,d)-motif which exists in each DNA.
 * hd: Hamming distance
 */
void search_motif(UINT32 hd)
{
    UINT32 j,k;
    UINT32 key,cmp_key;
    
    for(key=0;key<HASH_TABLE_LENGTH;key++){
        if(pKD_MOTIF[key]==NEW_PATTERN){
            pKD_MOTIF[key]=0;
            for(j=0;j<NUMBER_OF_DNAS;j++){
                for(k=0;k<=(SIZE_OF_DNA-KMER_LENGTH);k++){              
                    cmp_key = gen_hash_key(DNA_SRC+j*SIZE_OF_DNA+k,KMER_LENGTH);
                    if(hd>=hamming_distance_bykey(key,cmp_key,KMER_LENGTH)){
                        pKD_MOTIF[key]++;
                        //Check next DNA
                        break; 
                    }
                }
            }
        }
    }
    
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
    UINT32 i,j;
    UINT32 key;
    UINT32 match_cnt = 0;
    read_test_data("dataset_156_7.txt",DNA_SRC,MAX_ARRAY_SIZE);
#ifdef DEBUG
    printf("num:%d sz:%d kmer_sz:%d hd:%d \n",NUMBER_OF_DNAS,SIZE_OF_DNA,KMER_LENGTH,HAMMING_DISTANCE);
#endif
    // Init hash index.
    alpha['A'] = 0;
    alpha['C'] = 1;
    alpha['G'] = 2;
    alpha['T'] = 3;
    
    //Init hash table.
    pKD_MOTIF = malloc(HASH_TABLE_LENGTH);
    memset(pKD_MOTIF,0,HASH_TABLE_LENGTH*sizeof(UINT32));
    

    // Scan each DNA.
    for(i=0;i<NUMBER_OF_DNAS;i++){
        for(j=0;j<=(SIZE_OF_DNA-KMER_LENGTH);j++){
            
            key = gen_hash_key(DNA_SRC+j+i*SIZE_OF_DNA,KMER_LENGTH);
            //Add new key.
            if(0 == pKD_MOTIF[key])
                pKD_MOTIF[key] = NEW_PATTERN;
            //Generate key the hamming distance less than HAMMING_DISTANCE.          
            gen_kd_motif(key,0,HAMMING_DISTANCE,KMER_LENGTH);
            search_motif(HAMMING_DISTANCE);
        }
    }

    // Print out Motif
    for(i=0;i<HASH_TABLE_LENGTH;i++){
        if(pKD_MOTIF[i]==NUMBER_OF_DNAS){
            printpattern(i);
#ifdef DEBUG
            printf("%d\n",pKD_MOTIF[i]);
#endif            
        }
    }
}
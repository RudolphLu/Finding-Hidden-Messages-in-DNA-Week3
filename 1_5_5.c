/*
 * GreedyMotifSearch
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
#define SCAN_NUMOFDNA    1
#define SCAN_DNA         2



UINT8  DNA_SRC[MAX_ARRAY_SIZE];
UINT32 *PROFILE_PATTERN;
UINT32 pMostProbableKmer[2][100];//Store the most probability position.
UINT32 NUMBER_OF_DNA;
UINT32 SIZE_OF_DNA = 0;
UINT32 DNA_TOTAL_LENGTH;
UINT32 KMER_LENGTH = 0;
UINT32 alpha[256] = {0};


void read_test_data (const char* file_name,UINT8 *src ,UINT32 *profile,int maxsize) 
{    
    UINT32 phase = SCAN_KMER,count=0;
    UINT8 c;   
    FILE* file = fopen (file_name, "r");
    
    if (!file) return;

    NUMBER_OF_DNA = KMER_LENGTH = DNA_TOTAL_LENGTH = 0;
    
    do {
        if(SCAN_DNA == phase){
            c = fgetc(file);
            if( (c == 'G') || (c == 'C') || (c=='A') || (c=='T')){
                src[DNA_TOTAL_LENGTH ++] = c;                
            }else if((c == 0xd) || (c == 0xa)){
                continue;
            }
        }

        if(SCAN_KMER == phase){
            c = fgetc(file);
            if(c == ' '){
                phase = SCAN_NUMOFDNA;
                continue;
            }
            KMER_LENGTH *= 10;
            KMER_LENGTH += c-'0';
        }
        
        if(SCAN_NUMOFDNA == phase){
            c = fgetc(file);                     
            if((c == 0xd) || (c == 0xa)){
                if(0xa == c)
                    phase = SCAN_DNA;
                continue;
            }

            NUMBER_OF_DNA *= 10;
            NUMBER_OF_DNA += c-'0';
        }        

    } while (!feof (file) && count < maxsize);
    
    SIZE_OF_DNA = DNA_TOTAL_LENGTH/NUMBER_OF_DNA;
    
    fclose (file);
}

double caculate_pr(UINT8 *dna)
{
    UINT32 i;
    double multiple = 1;
    double j;
    
    for(i=0,multiple=1;i<KMER_LENGTH;i++){
        
        j = PROFILE_PATTERN[i]+PROFILE_PATTERN[i+KMER_LENGTH]+PROFILE_PATTERN[i+KMER_LENGTH*2]+PROFILE_PATTERN[i+KMER_LENGTH*3];
        j = PROFILE_PATTERN[alpha[dna[i]]*KMER_LENGTH+i]/j;
                    
        multiple*=j;
    }

    return multiple;
}




/*
 * Update MOTFIS probability table
 *    ----> Kmer
 * A  1 2
 * C  0 1
 * G  2 0
 * T  0 0
 */
void update_prtable(UINT8 *newmotfis)
{
    UINT32 i;
    for(i=0;i<KMER_LENGTH;i++){
        PROFILE_PATTERN[alpha[newmotfis[i]]*KMER_LENGTH+i]+=1;
    }
}

/*
 * Update MOTFIS probability table
 *    ----> Kmer
 * A  1 2
 * C  0 1
 * G  2 0
 * T  0 0
 */
void downgrade_prtable(UINT8 *newmotfis)
{
    UINT32 i;
    for(i=0;i<KMER_LENGTH;i++){
        PROFILE_PATTERN[alpha[newmotfis[i]]*KMER_LENGTH+i]--;
    }
}

void clean_prtable()
{
    UINT32 i;
    for(i=0;i<KMER_LENGTH;i++){
        PROFILE_PATTERN[i]=0;
        PROFILE_PATTERN[KMER_LENGTH+i]=0;
        PROFILE_PATTERN[2*KMER_LENGTH+i]=0;
        PROFILE_PATTERN[3*KMER_LENGTH+i]=0;
    }
}

void GreedyMotifSearch(UINT32 kmer,UINT32 t,UINT8 *dna_src)
{
    UINT32 score=0;
    UINT32 i=0,j,k;
    double m_prob,c_prob;
        
    for(j=1;j<t;j++){        
        for(i=0,m_prob=0,pMostProbableKmer[0][j] = 0;i<=(SIZE_OF_DNA-kmer);i++){              
              c_prob = caculate_pr(dna_src+i+j*SIZE_OF_DNA);
              if(c_prob>m_prob){
                  m_prob = c_prob;
                  //Update the max probabity position
                  pMostProbableKmer[0][j] = i;
              }              
        }
        k = pMostProbableKmer[0][j];
        //Update probability table.
        update_prtable(dna_src+j*SIZE_OF_DNA+k);
    }
}


/*
 * To get the minmum difference.
 */
void score(UINT8 *dna_src,UINT32 *best_score)
{
    UINT32 i,j,k,l;
    UINT32 sum_A,sum_C,sum_G,sum_T,max_NU=0,sum_column=0,sum=0;
    
    for(i=0;i<KMER_LENGTH;i++){
        
        sum_A = sum_C = sum_G = sum_T = 0;
        for(j=0;j<NUMBER_OF_DNA;j++){
            k=pMostProbableKmer[0][j];
            if(dna_src[j*SIZE_OF_DNA+k+i] == 'A') sum_A++;
            if(dna_src[j*SIZE_OF_DNA+k+i] == 'C') sum_C++;
            if(dna_src[j*SIZE_OF_DNA+k+i] == 'G') sum_G++;
            if(dna_src[j*SIZE_OF_DNA+k+i] == 'T') sum_T++;
        }
        
        sum_column = sum_A+sum_C+sum_G+sum_T;
        
        l = (sum_A>sum_C)?sum_A:sum_C;
        l = (sum_G>l)?sum_G:l;
        l = (sum_T>l)?sum_T:l;
        sum_column -= l;
        sum+=sum_column;
    }
    
    
    
    if(sum<*best_score){
        *best_score = sum;
        for(i=0;i<NUMBER_OF_DNA;i++){           
            pMostProbableKmer[1][i] = pMostProbableKmer[0][i];
        }
    }
    
    
}

void main()
{
    UINT32 i=0,j=0,k=0;
    UINT32 pr_multi,max_pr=0;
    UINT32 best_score=0xffffffff;
    read_test_data("dataset_158_9.txt",DNA_SRC,PROFILE_PATTERN,MAX_ARRAY_SIZE);
#ifdef DEBUG
    printf("num_d:%d dna_total:%d dna_sz:%d kmer:%d \n",NUMBER_OF_DNA,DNA_TOTAL_LENGTH,SIZE_OF_DNA,KMER_LENGTH);
#endif
    // Init hash index.
    alpha['A'] = 0;
    alpha['C'] = 1;
    alpha['G'] = 2;
    alpha['T'] = 3;
  
    //Init profile table
    PROFILE_PATTERN = malloc(sizeof(UINT32)*KMER_LENGTH*NUMBER_OF_DNA);
    memset(PROFILE_PATTERN,0,sizeof(UINT32)*KMER_LENGTH*NUMBER_OF_DNA);
    memset(pMostProbableKmer,0,sizeof(UINT32)*NUMBER_OF_DNA*2);


    for(i=0;i<=(SIZE_OF_DNA-KMER_LENGTH);i++){
        clean_prtable();
        update_prtable(DNA_SRC+i);        
        pMostProbableKmer[0][0]=i;
        GreedyMotifSearch(KMER_LENGTH,NUMBER_OF_DNA,DNA_SRC);
        // Score the min
        score(DNA_SRC,&best_score);  
    }

    for(i=0;i<NUMBER_OF_DNA;i++){
        j = pMostProbableKmer[1][i];
        k = i*SIZE_OF_DNA+j;
        for(j=0;j<KMER_LENGTH;j++)
            printf("%c",DNA_SRC[k+j]);
        
        printf("\n");
    }
  
}

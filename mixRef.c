/*
 * =====================================================================================
 *
 *       Filename:  mixRef.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/09/2012 01:25:07 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <errno.h>
#include "kseq.h"
#include "hapmap.h"
#include "mixRef.h"
#define CHAR_IN_WORD 8
/*
#if UNIT_TESTING
#include <stddef.h>
#include <setjmp.h>
#include <cmockery.h>
#endif 
*/
#include<stdint.h>
KSEQ_INIT(gzFile, gzread)
static uint8_t nt5_4bit_table[256] = {
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0 /*'-'*/, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 1, 0, 2,  0, 0, 0, 4,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  8, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 1, 0, 2,  0, 0, 0, 4,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  8, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
};// {A,C,G,T,#, N} = {0, 1, 2, 3, 4, 5}

#define CHAR_PER_BYTE 2
#define BYTE_PER_WORD 4
#define CHAR_PER_WORD 8
mixRef_t *mixRef_restore(const char *fn)
{
   
    FILE *fp;
    fp = fopen(fn,"rb");
    if(fp == NULL){
        fprintf(stderr, "[mixRef_restore]:file open %s fail!\n", fn);
        exit(EXIT_FAILURE);
    }

    int num_in_last_word; 
    mixRef_t *mixRef; 
    
    mixRef = calloc(1, sizeof(mixRef_t));
    if(mixRef == NULL){
        fprintf(stderr, "[mixRef]allocate mem fail!\n");
        exit(EXIT_FAILURE);
    } 
    fseek(fp, -4, SEEK_END);
    fread(&num_in_last_word, sizeof(uint32_t), 1, fp);
    mixRef->l =  (ftell(fp) -(2*BYTE_PER_WORD))*CHAR_PER_BYTE  + num_in_last_word;
    mixRef->seq = calloc((mixRef->l+CHAR_IN_WORD-1)/CHAR_IN_WORD, sizeof(uint32_t)); 
    if(mixRef->seq == NULL){
        fprintf(stderr, "[mixRef_restore]: mixRef->seq allocate mem fail!\n");
        exit(EXIT_FAILURE);
    }
    fseek(fp, 0, SEEK_SET);
    fread(mixRef->seq, sizeof(uint32_t),(mixRef->l + CHAR_IN_WORD -1)/CHAR_IN_WORD , fp); 
    fclose(fp);
    return mixRef;
}
void mixRef_destroy(mixRef_t *mixRef)
{
    free(mixRef->seq);
    free(mixRef);
}
#define __set_pac(pac, l, v) (pac)[(l)>>3]|=((v)<<4*((l)%8))
#define __get_pac(pac, l) (((pac)[(l)>>3]>>4*((l)%8))&15)
int build_mixRef(const char *fn_fa, const char *fn_hapmap, const char *fn_mixRef)
{
    int i;    
    FILE *fp_hapmap =NULL, *fp_mixRef = NULL;
    gzFile fp_fa = Z_NULL;
    kseq_t *seq; hapmap_t *hm;
    mixRef_t *mixRef; uint32_t mixRef_l;
    {//open stream
        fp_fa	= gzopen( fn_fa, "r" );
        if ( fp_fa == NULL ) {
            fprintf ( stderr, "couldn't open file '%s'; %s\n",
                fn_fa, strerror(errno) );
            exit (EXIT_FAILURE);
        }
        fp_hapmap	= fopen( fn_hapmap, "r" );
        if ( fp_hapmap == NULL ) {
            fprintf ( stderr, "couldn't open file '%s'; %s\n",
                    fn_hapmap, strerror(errno) );
            exit (EXIT_FAILURE);
        }
        fp_mixRef	= fopen( fn_mixRef, "w" );
        if ( fp_hapmap == NULL ) {
            fprintf ( stderr, "couldn't open file '%s'; %s\n",
                    fn_mixRef, strerror(errno) );
            exit (EXIT_FAILURE);
        }
    }
    
    mixRef = (mixRef_t *)calloc(1, sizeof(mixRef_t));
    if(mixRef == NULL){
        fprintf(stderr, "[main_gen_mixRef]allocate mem fail!\n");
        exit(EXIT_FAILURE);
    }
    seq = kseq_init(fp_fa);
    hm = hapmap_init(fp_hapmap); 
  
    uint32_t tot_l = 0; 
    int l = 0, l_last_word = 0;

    while( (l = kseq_read(seq)) >0)
    {
        int i_seq = 0;

        mixRef_l = mixRef->l/CHAR_IN_WORD *CHAR_IN_WORD +1;
        if(tot_l +l > mixRef_l){
//            mixRef->seq =   realloc(mixRef->seq, (tot_l+l)/CHAR_IN_WORD *CHAR_IN_WORD +1);
            uint32_t seq_size_cur_in_byte, seq_size_last_in_byte, seq_size_cur_in_word,seq_size_last_in_word; 
            seq_size_cur_in_word = (tot_l +l +CHAR_IN_WORD -1)/CHAR_IN_WORD;
            seq_size_last_in_word = (tot_l + CHAR_IN_WORD -1)/CHAR_IN_WORD;

            seq_size_cur_in_byte = seq_size_cur_in_word*BYTE_PER_WORD;     
            seq_size_last_in_byte = seq_size_last_in_word*BYTE_PER_WORD;
            mixRef->seq =   realloc(mixRef->seq, seq_size_cur_in_byte);     
//            memset(mixRef->seq + mixRef->l/CHAR_IN_WORD*CHAR_IN_WORD+1, (tot_l+l)/CHAR_IN_WORD *CHAR_IN_WORD - mixRef->l/CHAR_IN_WORD*CHAR_IN_WORD, 0); 
//            memset(mixRef->seq);
            memset(mixRef->seq + seq_size_last_in_word, seq_size_cur_in_byte - seq_size_last_in_byte, 0);
            if(mixRef->seq  == NULL){
                fprintf(stderr, "[main_gen_mixRef]:realloc mem fail!\n");
                exit(EXIT_FAILURE);
            }
            mixRef->l = tot_l +l;
        } 
/****************************************************************/        
        uint32_t *start, *end;
        uint32_t word;
        
        start = mixRef->seq + tot_l/CHAR_IN_WORD;            
        end = mixRef->seq + mixRef->l/CHAR_IN_WORD;
        //first word
        /*
        word = 0;
        for(i = 0; i < l_last_word; ++i){
            word |= (*start>>i*4) &15;
        }*/
        if(tot_l > 0){
            word = *start;
            for(i = l_last_word; i < CHAR_IN_WORD; ++i){
                word|= nt5_4bit_table[seq->seq.s[i_seq]]<< 4*i; 
                ++i_seq;
            }
            *start = word;
            ++start;
        }        
        //core loop
        for(;start!= end; ++start)
        {
            word = 0;
            for(i = 0; i < CHAR_IN_WORD; ++i){
                word |= nt5_4bit_table[seq->seq.s[i_seq]]<<4*i; 
                ++i_seq;
#ifdef DEBUG
fprintf(stderr, "%x\n", word);
#endif
            }
           *start = word; 
        }
        //last word
        l_last_word = mixRef->l % CHAR_IN_WORD;
        word = 0;
        for(i = 0; i < l_last_word; ++i){
            word |= nt5_4bit_table[seq->seq.s[i_seq]]<<4*i; 
#ifdef DEBUG
fprintf(stderr, "%x\n", word);
#endif
            ++i_seq;
        }
        *start = word;
/***************************************************************/       
        hapmap_readhm(hm);
fprintf(stderr, "snppos:%u\n", hm->snp_pos[0]);
        for(i=0; i < hm->snp_num; ++i){
             __set_pac(mixRef->seq, tot_l+hm->snp_pos[i], hm->snp_type[i] &15); 
        }        
/***************************************************************/        
        tot_l += l;    
    } 
#ifdef DEBUG
int k;
for(k = 0; k < mixRef->l/CHAR_IN_WORD +1; ++k){
    fprintf(stderr, "%x\n", mixRef->seq[k]);
}
#endif
    fwrite(mixRef->seq, sizeof(uint32_t), (mixRef->l+1-CHAR_IN_WORD)/CHAR_IN_WORD, fp_mixRef);
    
    uint32_t count_in_last_word = mixRef->l %CHAR_PER_WORD;
    if(count_in_last_word == 0){
        fwrite(&count_in_last_word, sizeof(uint32_t), 1, fp_mixRef);
    }
    fwrite(&count_in_last_word, sizeof(uint32_t), 1, fp_mixRef);
    
    free(mixRef->seq);
    free(mixRef);
    kseq_destroy(seq);
    hapmap_destroy(hm); 
    
    
    {//close stream
        if( gzclose(fp_fa) == EOF ) {			/* close output file   */
            fprintf ( stderr, "couldn't close file '%s'; %s\n",
                fn_fa, strerror(errno) );
            exit (EXIT_FAILURE);
        }
        if( fclose(fp_hapmap) == EOF ) {			/* close input file   */
            fprintf ( stderr, "couldn't close file '%s'; %s\n",
                    fn_hapmap, strerror(errno) );
            exit (EXIT_FAILURE);
        }
        if( fclose(fp_mixRef) == EOF ) {			/* close input file   */
            fprintf ( stderr, "couldn't close file '%s'; %s\n",
                    fn_hapmap, strerror(errno) );
            exit (EXIT_FAILURE);
        }
    }
    return EXIT_SUCCESS;
}
#ifdef MAIN_GEN_MIXREF
static int usage()
{
    fprintf(stderr, "********************\n");
    fprintf(stderr, "gm <ref.fa.in> <hapmap.in><mixRef.out>\n");
     
    fprintf(stderr, "********************\n");
    return 1;
}

int main(int argc, char *argv[])
{
    if(argc < 4){
        return usage();
    }
    
    char *fn_fa = argv[1];
    char *fn_hapmap = argv[2];
    char *fn_mixRef = argv[3];
    return build_mixRef(fn_fa, fn_hapmap, fn_mixRef); 
}
#endif 

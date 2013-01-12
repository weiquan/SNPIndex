/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Quan Wei <wquanhit@gmail.com> */
/*改为用4bit表示每个字符*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "bntseq.h"

#include "utils.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

unsigned char nst_nt5_table[256] = {
	7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7, 
	7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7, 
	7, 7, 7, 4,  7, 7, 7, 7,  7, 7, 7, 7,  7, 7 /*'-'*/, 7, 7,
	7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7, 
	7, 0, 7, 1,  7, 7, 7, 2,  7, 7, 7, 7,  7, 7, 5, 7, 
	7, 7, 7, 7,  3, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7, 
	7, 0, 7, 1,  7, 7, 7, 2,  7, 7, 7, 7,  7, 7, 5, 7, 
	7, 7, 7, 7,  3, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7, 
	7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7, 
	7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7, 
	7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7, 
	7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7, 
	7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7, 
	7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7, 
	7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7, 
	7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7,  7, 7, 7, 7
};// {A,C,G,T,#, N} = {0, 1, 2, 3, 4, 5}

void R_bns_dump(const bntseq_t *bns, const char *prefix)
{
	char str[1024];
	FILE *fp;
	int i;
	{ // dump .ann
		strcpy(str, prefix); strcat(str, ".ann");
		fp = xopen(str, "w");
		fprintf(fp, "%lld %d %u\n", (long long)bns->l_pac, bns->n_seqs, bns->seed);
		for (i = 0; i != bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			fprintf(fp, "%d %s", p->gi, p->name);
			if (p->anno[0]) fprintf(fp, " %s\n", p->anno);
			else fprintf(fp, "\n");
			fprintf(fp, "%lld %d %d\n", (long long)p->offset, p->len, p->n_ambs);
		}
		fclose(fp);
	}
	{ // dump .amb
		strcpy(str, prefix); strcat(str, ".amb");
		fp = xopen(str, "w");
		fprintf(fp, "%lld %d %u\n", (long long)bns->l_pac, bns->n_seqs, bns->n_holes);
		for (i = 0; i != bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			fprintf(fp, "%lld %d %c\n", (long long)p->offset, p->len, p->amb);
		}
		fclose(fp);
	}
}

bntseq_t *R_bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename)
{
	char str[1024];
	FILE *fp;
	bntseq_t *bns;
	long long xx;
	int i;
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	{ // read .ann
		fp = xopen(ann_filename, "r");
		fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
		bns->l_pac = xx;
		bns->anns = (bntann1_t*)calloc(bns->n_seqs, sizeof(bntann1_t));
		for (i = 0; i < bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			char *q = str;
			int c;
			// read gi and sequence name
			fscanf(fp, "%u%s", &p->gi, str);
			p->name = strdup(str);
			// read fasta comments 
			while ((c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;
			*q = 0;
			if (q - str > 1) p->anno = strdup(str + 1); // skip leading space
			else p->anno = strdup("");
			// read the rest
			fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
			p->offset = xx;
		}
		fclose(fp);
	}
	{ // read .amb
		int64_t l_pac;
		int32_t n_seqs;
		fp = xopen(amb_filename, "r");
		fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
		l_pac = xx;
		xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
		bns->ambs = (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t));
		for (i = 0; i < bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			fscanf(fp, "%lld%d%s", &xx, &p->len, str);
			p->offset = xx;
			p->amb = str[0];
		}
		fclose(fp);
	}
	{ // open .pac
		bns->fp_pac = xopen(pac_filename, "rb");
	}
	return bns;
}

bntseq_t *R_bns_restore(const char *prefix)
{  
	char ann_filename[1024], amb_filename[1024], pac_filename[1024];
	strcat(strcpy(ann_filename, prefix), ".ann");
	strcat(strcpy(amb_filename, prefix), ".amb");
	strcat(strcpy(pac_filename, prefix), ".pac");
	return R_bns_restore_core(ann_filename, amb_filename, pac_filename);
}

void R_bns_destroy(bntseq_t *bns)
{
	if (bns == 0) return;
	else {
		int i;
		if (bns->fp_pac) fclose(bns->fp_pac);
		free(bns->ambs);
		for (i = 0; i < bns->n_seqs; ++i) {
			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		free(bns->anns);
		free(bns);
	}
}
#define MASK_4BIT 15
#define N_nt 5
#define NT_PER_BYTE 2
#define BITS_PER_NT 4
#define _set_pac(pac, l, c) ((pac)[(l)>>1] |= (c)<<((~(l)&1)<<2))
#define _get_pac(pac, l) ((pac)[(l)>>1]>>((~(l)&1)<<2)&MASK_4BIT)

static uint8_t *R_add1(const kseq_t *seq, bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q)
{
	bntann1_t *p;
	int i, lasts;
	if (bns->n_seqs == *m_seqs) {
		*m_seqs <<= 1;
		bns->anns = (bntann1_t*)realloc(bns->anns, *m_seqs * sizeof(bntann1_t));
	}
	p = bns->anns + bns->n_seqs;
	p->name = strdup((char*)seq->name.s);
	p->anno = seq->comment.s? strdup((char*)seq->comment.s) : strdup("(null)");
	p->gi = 0; p->len = seq->seq.l;
	p->offset = (bns->n_seqs == 0)? 0 : (p-1)->offset + (p-1)->len;
	p->n_ambs = 0;
	for (i = lasts = 0; i < seq->seq.l; ++i) {
		int c = nst_nt5_table[(int)seq->seq.s[i]];
		if (c >= N_nt) { // N
			if (lasts == seq->seq.s[i]) { // contiguous N
				++(*q)->len;
			} else {
				if (bns->n_holes == *m_holes) {
					(*m_holes) <<= 1;
					bns->ambs = (bntamb1_t*)realloc(bns->ambs, (*m_holes) * sizeof(bntamb1_t));
				}
				*q = bns->ambs + bns->n_holes;
				(*q)->len = 1;
				(*q)->offset = p->offset + i;
				(*q)->amb = seq->seq.s[i];
				++p->n_ambs;
				++bns->n_holes;
			}
		}
		lasts = seq->seq.s[i];
		{ // fill buffer
			if (c >= N_nt) c = lrand48()&3;
			if (bns->l_pac == *m_pac) { // double the pac size
				*m_pac <<= 1;
				pac = realloc(pac, *m_pac/NT_PER_BYTE);
				memset(pac + bns->l_pac/NT_PER_BYTE, 0, (*m_pac - bns->l_pac)/NT_PER_BYTE);
			}
			_set_pac(pac, bns->l_pac, c);
			++bns->l_pac;
		}
	}
	++bns->n_seqs;
	return pac;
}

int64_t R_bns_fasta2bntseq(gzFile fp_fa, const char *prefix)
{
//	extern void seq_reverse(int len, ubyte_t *seq, int is_comp); // in bwaseqio.c
	kseq_t *seq;
	char name[1024];
	bntseq_t *bns;
	uint8_t *pac = 0;
    uint8_t *reverse_pac = NULL;
	int32_t m_seqs, m_holes;
	int64_t ret = -1, m_pac, l, m_r_pac;
	bntamb1_t *q;
	FILE *fp, *fp_r;
    int i;
	// initialization
	seq = kseq_init(fp_fa);
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	bns->seed = 11; // fixed seed for random generator
	srand48(bns->seed);
	m_seqs = m_holes = 8; m_pac = 0x10000;
	bns->anns = (bntann1_t*)calloc(m_seqs, sizeof(bntann1_t));
	bns->ambs = (bntamb1_t*)calloc(m_holes, sizeof(bntamb1_t));
	pac = calloc(m_pac/NT_PER_BYTE, 1);
	q = bns->ambs;
	strcpy(name, prefix); strcat(name, ".pac");
	fp = xopen(name, "wb");
    memset(name, '\0', 1024);strcpy(name, prefix); strcat(name, ".rpac");
    fp_r = xopen(name, "wb");
	// read sequences
	while (kseq_read(seq) >= 0) pac = R_add1(seq, bns, pac, &m_pac, &m_seqs, &m_holes, &q);
    /*
    if (!for_only) { // add the reverse complemented sequence
		m_pac = (bns->l_pac * 2 + 3) / 4 * 4;
		pac = realloc(pac, m_pac/4);
		memset(pac + (bns->l_pac+3)/4, 0, (m_pac - (bns->l_pac+3)/4*4) / 4);
		for (l = bns->l_pac - 1; l >= 0; --l, ++bns->l_pac)
			_set_pac(pac, bns->l_pac, 3-_get_pac(pac, l));
	}
	ret = bns->l_pac;*/
    fprintf(stderr, "[bns_fasta2bntseq]:reverse_pac!\n");
    m_r_pac = (bns->l_pac +NT_PER_BYTE-1)/NT_PER_BYTE *NT_PER_BYTE;
    reverse_pac = calloc(m_r_pac/NT_PER_BYTE, sizeof(uint8_t));
    for(l = bns->l_pac-1, i =0; l>=0, i < bns->l_pac; --l, ++i)
    {
        _set_pac(reverse_pac, i, _get_pac(pac, l)); 
    } 
	ubyte_t ct;
    { // finalize .pac and  file
		fwrite(pac, 1, (bns->l_pac>>1) + ((bns->l_pac&1) == 0? 0 : 1), fp);
		// the following codes make the pac file size always (l_pac/4+1+1)
		if (bns->l_pac % NT_PER_BYTE == 0) {
			ct = 0;
			fwrite(&ct, 1, 1, fp);
		}
		ct = bns->l_pac % NT_PER_BYTE;
		fwrite(&ct, 1, 1, fp);
		// close .pac file
		fclose(fp);
	}	
    { // finalize .rpac and  file
		fwrite(reverse_pac, 1, (bns->l_pac>>1) + ((bns->l_pac&1) == 0? 0 : 1), fp_r);
		// the following codes make the pac file size always (l_pac/4+1+1)
		if (bns->l_pac % NT_PER_BYTE == 0) {
			ct = 0;
			fwrite(&ct, 1, 1, fp_r);
		}
		ct = bns->l_pac % NT_PER_BYTE;
		fwrite(&ct, 1, 1, fp_r);
		// close .pac file
		fclose(fp_r);
	}
	bns_dump(bns, prefix);
	bns_destroy(bns);
	kseq_destroy(seq);
	free(pac);
    free(reverse_pac);
	return ret;
}

int R_fa2pac(int argc, char *argv[])
{
	int c, for_only = 0;
	gzFile fp;
	while ((c = getopt(argc, argv, "f")) >= 0) {
		switch (c) {
			case 'f': for_only = 1; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: fa2pac [-f] <in.fasta> [<out.prefix>]\n");
		return 1;
	}
	fp = xzopen(argv[optind], "r");
	R_bns_fasta2bntseq(fp, (optind+1 < argc)? argv[optind+1] : argv[optind]);
	gzclose(fp);
	return 0;
}
/*
int main (int argc, char *argv[])
{

    return fa2pac(argc, argv);
}*/
int R_bns_cnt_ambi(const bntseq_t *bns, int64_t pos_f, int len, int *ref_id)
{
	int left, mid, right, nn;
	if (ref_id) {
		left = 0; mid = 0; right = bns->n_seqs;
		while (left < right) {
			mid = (left + right) >> 1;
			if (pos_f >= bns->anns[mid].offset) {
				if (mid == bns->n_seqs - 1) break;
				if (pos_f < bns->anns[mid+1].offset) break; // bracketed
				left = mid + 1;
			} else right = mid;
		}
		*ref_id = mid;
	}
	left = 0; right = bns->n_holes; nn = 0;
	while (left < right) {
		mid = (left + right) >> 1;
		if (pos_f >= bns->ambs[mid].offset + bns->ambs[mid].len) left = mid + 1;
		else if (pos_f + len <= bns->ambs[mid].offset) right = mid;
		else { // overlap
			if (pos_f >= bns->ambs[mid].offset) {
				nn += bns->ambs[mid].offset + bns->ambs[mid].len < pos_f + len?
					bns->ambs[mid].offset + bns->ambs[mid].len - pos_f : len;
			} else {
				nn += bns->ambs[mid].offset + bns->ambs[mid].len < pos_f + len?
					bns->ambs[mid].len : len - (bns->ambs[mid].offset - pos_f);
			}
			break;
		}
	}
	return nn;
}
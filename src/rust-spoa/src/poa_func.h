#ifndef POA_FUNC_H
#define POA_FUNC_H

#ifdef __cplusplus
extern "C" {
#endif
unsigned poa_func(char** seqs,        // the sequences (null-terminated) to perform multiple-sequence-alignment with.
                  char** quals,       // phred quality scores.
                  int num_seqs,       // the number of sequences being multiply aligned
                  char* consensus,    // this chunk of memory will hold the return value, the consensus of the multiple alignment
                  int consensus_len,  // the amount of memory allocated to consensus, i.e. use the MAXIMUM length of the consensus
                  int alignment_type, // the alignment type: 0 = local align, 1 = global align, 2 = semi-global
                  int match_score,    // the score to give a sequence match in alignment, e.g. 5
                  int mismatch_score, // the score to give a sequence mismatch in alignment, e.g. -4
                  int gap_open,
                  int gap_extend);     // the score to give a sequence gap in alignment, e.g. -8

#ifdef __cplusplus
}
#endif

#endif // POA_FUNC_H

#include "poa_func.h"
#include "spoa/spoa.hpp"

extern "C" {

    // see the C header file (poa_func.h) for detailed descriptions of each argument
    unsigned poa_func(char** seqs, char** quals, int num_seqs,
                      char* consensus, int consensus_len,
                      int alignment_type, int match_score, int mismatch_score, int gap_open, int gap_extend) {

        if (num_seqs == 0) {
            return (unsigned) 0;
        }

        // populate the list of sequences
        std::vector<std::string> sequences;
        std::vector<std::string> quality_strings;
        size_t max_seq_len = 0;
        for (int i = 0; i < num_seqs; i++){
            sequences.push_back((std::string) seqs[i]);
            quality_strings.push_back((std::string) quals[i]);
            max_seq_len = std::max(max_seq_len, sequences.back().length());
        }

        auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(alignment_type),
                                                            (int8_t) match_score,
                                                            (int8_t) mismatch_score,
                                                            (int8_t) gap_open,
                                                            (int8_t) gap_extend);
        auto graph = spoa::createGraph();

        // add each of the real sequences (e.g. noisy sequence reads) to the graph

        for (int i = 0; i < num_seqs; i++) {
            //GOOD
            //auto alignment = (*alignment_engine)(sequences[i], graph);

            //TESTING
            auto alignment = alignment_engine->align(sequences[i], graph);
            graph->add_alignment(alignment, sequences[i], quality_strings[i]);
        }

        // generate the consensus sequence, assign it to the allocated memory block, and return the consensus length.
        std::string cns = graph->generate_consensus();

        int l = cns.length();
        if (l > consensus_len) {
            l = consensus_len;
        }

        for (int i = 0; i < l; i++){
            consensus[i] = cns[i];
        }

        sequences.clear();

        return (unsigned) l;
    }
}

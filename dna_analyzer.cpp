#include "dna_analyzer.h"
#include <sstream>
#include <iomanip>
#include <cmath>

// Embedding generator implementation
SequenceEmbedding get_sequence_embedding(const std::map<KmerKey,int>& kmer_counts, const std::string& sequence) {
    SequenceEmbedding emb;
    int total = sequence.length();
    int gc_count = 0, at_count = 0, rare_count = 0, pal_count = 0;
    double entropy = 0.0;

    // GC and AT content
    for (char c : sequence) {
        if (c == 'G' || c == 'C') gc_count++;
        else if (c == 'A' || c == 'T') at_count++;
    }

    // Rare k-mers, entropy, palindromes (for k=5)
    for (const auto& p : kmer_counts) {
        int count = p.second;
        std::string kmer = ""; // You need a way to reconstruct string from KmerKey if possible
        if (count < 2) rare_count += count;
        // For simplicity, use entropy = -sum(p*log2(p)), here p = count/total_kmers
        double freq = double(count)/kmer_counts.size();
        if (freq > 0) entropy -= freq * log2(freq);
        // If you can reconstruct kmer: if (kmer == reverse(kmer)) pal_count += count;
    }

    int total_kmers = 0;
    for (auto& p : kmer_counts) total_kmers += p.second;

    emb.vector[0] = double(gc_count)/total;          // GC content
    emb.vector[1] = double(rare_count)/total_kmers;  // Fraction of rare k-mers
    emb.vector[2] = double(at_count)/total;          // AT content
    emb.vector[3] = entropy/kmer_counts.size();      // Normalized entropy
    emb.vector[4] = double(pal_count)/total_kmers;   // Palindromic k-mers

    // Normalize embedding
    double mag = 0.0;
    for (int i=0;i<emb.size;i++) mag += emb.vector[i]*emb.vector[i];
    mag = std::sqrt(mag);
    if (mag > 0.0) for (int i=0;i<emb.size;i++) emb.vector[i] /= mag;

    return emb;
}


// Prediction implementation
std::string predict_function_from_embedding(const SequenceEmbedding& emb, double gc) {
    const double R[5] = {0.9, 0.1, 0.1, 0.9, 0.2}; // Ribosomal: High GC, high entropy
    const double T[5] = {0.1, 0.9, 0.1, 0.1, 0.9};
    // Metabolic: balanced GC/AT, moderate entropy
    const double M[5] = {0.5, 0.2, 0.6, 0.7, 0.1};

    double r = 0.0, t = 0.0, m = 0.0;
    for (int i = 0; i < emb.size; ++i) {
        r += emb.vector[i] * R[i];
        t += emb.vector[i] * T[i];
        m += emb.vector[i] * M[i];
    }
    double total_score = r + t + m;
    double r_frac = r / total_score;
    double t_frac = t / total_score;
    double m_frac = m / total_score;

    std::ostringstream ss;
    ss << std::fixed << std::setprecision(4);
    ss << "Scores: [Ribosomal: " << r << ", Mobile: " << t << ", Metabolic: " << m << "]\n";
    ss << "Fractions: [Ribosomal: " << r_frac << ", Mobile: " << t_frac << ", Metabolic: " << m_frac << "]\n";

    if (r_frac > 0.60) {
        ss << "PREDICTION: Ribosomal / Structural (dominant)";
    } else if (t_frac > 0.60) {
        ss << "PREDICTION: Mobile Element / Repeats (dominant)";
    } else if (m_frac > 0.60) {
        ss << "PREDICTION: Metabolic / Housekeeping (dominant)";
    } else if (r_frac > 0.40 && r_frac > t_frac && r_frac > m_frac) {
        ss << "PREDICTION: Ribosomal (mixed signatures)";
    } else if (t_frac > 0.40 && t_frac > r_frac && t_frac > m_frac) {
        ss << "PREDICTION: Mobile Element (mixed signatures)";
    } else if (m_frac > 0.40 && m_frac > r_frac && m_frac > t_frac) {
        ss << "PREDICTION: Metabolic (mixed signatures)";
    } else {
        ss << "PREDICTION: No single function dominant, potential mosaic or unknown function.";
    }

    return ss.str();
}


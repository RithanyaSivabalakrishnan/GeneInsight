#ifndef DNA_ANALYZER_H
#define DNA_ANALYZER_H

#include <vector>
#include <string>
#include <map>
#include <utility>
#include <stdexcept>
#include <limits>
#include <cmath>

using KmerKey = unsigned long long;

struct SequenceEmbedding {
    double vector[5];
    const int size = 5;
    SequenceEmbedding() { for (int i=0;i<size;++i) vector[i]=0.0; }
};

// KMerHashTable definition (kept in header - class only)
class KMerHashTable {
private:
    const KmerKey EMPTY_KEY = std::numeric_limits<KmerKey>::max();
    std::vector<std::pair<KmerKey, std::vector<int>>> table;
    int size_;
    int K;
    int unique_keys_count;
    const double MAX_OCCUPANCY = 0.7;

    inline KmerKey base_to_int(char base) const {
        switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: throw std::runtime_error("Invalid base in k-mer");
        }
    }

    KmerKey encode_kmer(const std::string &kmer) const {
        if ((int)kmer.size() != K) throw std::runtime_error("k length mismatch");
        KmerKey key = 0;
        for (char c : kmer) {
            key = (key << 2) | base_to_int(c);
        }
        return key;
    }

    int hash(KmerKey key) const { return (int)(key % (KmerKey)size_); }

    int find_slot_for_insert(KmerKey key) const {
        int idx = hash(key);
        int start = idx;
        while (table[idx].first != EMPTY_KEY && table[idx].first != key) {
            idx = (idx + 1) % size_;
            if (idx == start) return -1;
        }
        return idx;
    }

    void rehash(int new_size) {
        std::vector<std::pair<KmerKey, std::vector<int>>> old = std::move(table);
        int old_size = size_;
        table.assign(new_size, {EMPTY_KEY, {}});
        size_ = new_size;
        unique_keys_count = 0;
        for (int i=0;i<old_size;++i) {
            if (old[i].first != EMPTY_KEY) {
                int slot = find_slot_for_insert(old[i].first);
                if (slot == -1) throw std::runtime_error("rehash failure");
                table[slot].first = old[i].first;
                table[slot].second = std::move(old[i].second);
                unique_keys_count++;
            }
        }
    }

    void insert_internal(KmerKey key, int pos) {
        if ((double)unique_keys_count / (double)size_ > MAX_OCCUPANCY) {
            rehash(size_ * 2);
        }
        int slot = find_slot_for_insert(key);
        if (slot == -1) throw std::runtime_error("Hash table full");
        if (table[slot].first == EMPTY_KEY) {
            table[slot].first = key;
            unique_keys_count++;
        }
        table[slot].second.push_back(pos);
    }

public:
    KMerHashTable(int initial_size, int k) : size_(initial_size), K(k), unique_keys_count(0) {
        if (K <= 0 || K > 32) throw std::runtime_error("K must be between 1 and 32");
        table.assign(size_, {EMPTY_KEY, {}});
    }

    void build_index(const std::string &genome) {
        table.assign(size_, {EMPTY_KEY, {}});
        unique_keys_count = 0;
        if ((int)genome.size() < K) return;
        for (int i=0;i<= (int)genome.size()-K; ++i) {
            std::string kmer = genome.substr(i, K);
            try {
                KmerKey key = encode_kmer(kmer);
                insert_internal(key, i);
            } catch (...) { /* skip invalid */ }
        }
    }

    std::vector<int> query(const std::string &kmer_query) const {
        if ((int)kmer_query.size() != K) return {};
        try {
            KmerKey key = encode_kmer(kmer_query);
            int idx = hash(key);
            int start = idx;
            while (table[idx].first != EMPTY_KEY) {
                if (table[idx].first == key) return table[idx].second;
                idx = (idx + 1) % size_;
                if (idx == start) break;
            }
        } catch (...) {}
        return {};
    }

    std::map<KmerKey,int> get_kmer_counts() const {
        std::map<KmerKey,int> counts;
        for (const auto &p : table) {
            if (p.first != EMPTY_KEY) counts[p.first] = (int)p.second.size();
        }
        return counts;
    }

    int get_K() const { return K; }
    int get_unique_count() const { return unique_keys_count; }
};

// Declare embedding & prediction functions (defined in .cpp)
SequenceEmbedding get_sequence_embedding(const std::map<KmerKey, int>& kmer_counts, const std::string& sequence);
std::string predict_function_from_embedding(const SequenceEmbedding &embedding, double gc_content);

#endif // DNA_ANALYZER_H


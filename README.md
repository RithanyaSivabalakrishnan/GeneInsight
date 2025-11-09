# GeneInsight
A simple, interactive DNA sequence toolkit for k-mer based search, structure analysis, and genome function prediction using language-inspired embeddings and fast hash tables.

# Overview
KmerCanvas is designed to help you:

~ Search for any k-mer in large DNA sequences and view all their positions instantly.

~ Predict whether a sequence is ribosomal, metabolic, or a repetitive/mobile element using a compact embedding and scoring approach inspired by natural language processing (NLP).

# Features
~ K-mer Search: Rapidly find and count any k-mer (substring of length k) in a genomic sequence.

~ Embedding-Based Prediction: Converts sequences into numerical vectors capturing GC/AT content, rarity, entropy, and more—then predicts the likely biological function.

~ GC Content: Computes percentage GC content for any sequence or region.

~ Intuitive GUI: Open and analyze FASTA files with a clean, user-friendly interface.

# NLP Concepts Applied
~ Tokenization: Splits DNA into k-mers (analogous to words or n-grams).

~ Feature Extraction: Generates feature vectors for each sequence (GC content, AT content, sequence entropy, rare and palindromic k-mers).

~ Embedding: Converts genomic “language” to a vector, similar to text embeddings in NLP tasks.

~ Classification: Maps embedding to function classes (ribosomal, mobile/repeat, metabolic) using weighted scoring.

# Hashmap Concepts
~ K-mer Hashing: Each k-mer is encoded as a unique integer using a custom base-4 converter and stored in a hash table (fast lookup, low memory).

~ Open Addressing: Resolves hash collisions, ensuring efficient use of space and guaranteeing quick access.

~ Position Mapping & Counting: Tracks k-mer occurrence and all positions in sequence, supporting advanced queries and further embedding computation.


# Data Credits
Sequence data used and tested with this tool were obtained from the NCBI GenBank/RefSeq/ENA database(s).

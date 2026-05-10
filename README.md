# DNA Motif Finder: A Graph-Based Approach

### An algorithm for discovering common DNA motifs in biological sequences, integrating read quality analysis (Phred scores).

## About the Project
This C++ program identifies recurring substrings (motifs) across multiple DNA sequences simultaneously. The project leverages graph theory, representing potential motifs as vertices and their similarities as edges. The primary objective is to find a **clique**, which represents a specific motif present in every analyzed sequence.

## Key Features
- **Quality Control (QC) Filter:** The program processes `.qual` file and discards nucleotides with low reliability scores (below a user-defined threshold).
- **Nucleotide Positioning:** Maintains original sequence indexing even after filtering out low-quality bases, ensuring biological accuracy.
- **Graph Optimization (Hash Maps):** Graph construction utilizes (`unordered_map`) to drastically accelerate the connection process between identical oligonucleotides.
- **Clique Detection:** The algorithm identifies a set of interconnected vertices where each vertex originates from a different input sequence.
- **Biological Validation:** Accounts for the physical distance between nucleotides using a proximity parameter (`10 * dl_podciag`).

## Technologies
- **Language:** C++
- **Data Structures:** Graphs (adjacency lists), Vectors, Associative Maps (Hash Maps).
- **Algorithms:** Degree-based vertex sorting, graph clique detection.

## Input File Structure
The program requires two files in the same directory:
1. `sekwencje.fasta` – DNA sequences in FASTA format (example available)
2. `sekwencje.qual` – Quality values (Phred scores) for each base (example available)

## Compilation and Execution

### Compilation:
Use the `g++` compiler (the `-O3` optimization flag is highly recommended for large datasets):
```bash
g++ -O3 main.cpp -o motif_finder
```
### Execution:
```bash
./motif_finder

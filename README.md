# Klebsiella phage RBP Atlas

An interactive browser-based atlas for exploring the structural and functional diversity of **receptor-binding proteins (RBPs)** from *Klebsiella* bacteriophages. The atlas integrates predicted protein structures, ECOD domain annotations, pseudo-domain classifications, and sequence modularity evidence into a unified interactive visualization.

---

## Overview

RBPs are the host-recognition modules of bacteriophages, determining which bacterial capsule serotypes a phage can infect. This atlas presents **382 RBPs** from *Klebsiella* phages, classified into **39 structural classes** and organized into sequence and structural clusters. It covers:

- Domain architecture annotation via ECOD (Evolutionary Classification of Protein Domains)
- 3D structure rendering with domain-matched coloring
- **100 pseudo-domain clusters** mapping structurally conserved regions across RBP classes
- Sequence modularity evidence across **24 RBP classes** from **2,489 protein-pair alignments**

---

## Features

### Six explorer views

| View | Description |
|---|---|
| **RBPs** | All 382 proteins grouped by phage morphology |
| **RBP Classes** | Proteins grouped by structural class (39 classes) |
| **RBP Clusters** | Proteins grouped by sequence cluster |
| **ECOD Domains** | Three-level hierarchy: domain type → class → proteins, with domain bar rows |
| **Pseudo-domain Clusters** | 100 structural pseudo-domain clusters mapped onto RBP structures with TM-score badges |
| **Sequence Modularity** | Interactive Cytoscape.js network: nodes = RBP classes (size ∝ membership), edges = modularity evidence colored by category |

### Detail panel
- **Interactive 3D structure** — NGL-powered viewer with B-factor, ECOD, or pseudo-domain coloring; spin, reset, surface overlay, and fullscreen mode
- **Sequence coverage strip** — Canvas-drawn strip with hover-to-highlight: hovering over a residue position highlights the corresponding atom in the 3D viewer
- **Domain architecture** — Large domain bar with legend (ECOD view)
- **Pseudo-domain region** — Proportional region bar with TM-score metadata (PD view)
- **Protein & genome metadata** — Functional annotations (Pfam, ECOD, PHROGS), prediction scores, phage taxonomy, genome size, capsule type, host range, and PubMed links
- **Download** — `.pdb.gz` structure file download button

### Sequence modularity network
- Nodes sized proportional to class membership
- Edges colored by modularity category:
  - **Cterminal variation** — blue (`#1f77b4`)
  - **Nterminal sharing** — orange (`#ff7f0e`)
  - **Putative recombination hotspots** — green (`#2ca02c`)
  - **Cterminal sharing** — red (`#d62728`)
  - **Conserved regions** — purple (`#9467bd`)
- Click a **node** → highlight neighborhood, show protein-pair interaction table for that class
- Click an **edge** → highlight that class-pair connection, show protein pairs for that specific category
- Export TSV: all interactions, per-node, or per-edge

### Other
- **TSV export** — Download the current view's data as a tab-separated file
- **Universal search** — Searches across protein ID, phage name, genome name, genus, capsule type (KL1/KL2/…), phage cluster, morphology, host range, RBP class, RBP cluster, pseudo-domain cluster, and ECOD domain; results grouped by category with direct navigation
- **Resizable detail panel** — Drag the left edge to resize; panel supports vertical resize of the NGL viewport
- **Home donut chart** — ECOD domain coverage across all structured proteins

---

## Data files

| File | Description |
|---|---|
| `rbp_table_v16122025.tsv` | RBP metadata: functional predictions, model quality, cluster/class assignments |
| `ecod-map.tsv` | ECOD domain hits mapped to each RBP |
| `genome_table_v25122025.tsv` | Phage genome metadata: taxonomy, morphology, host range, capsule, publications |
| `ecod-color-map.txt` | Color assignments for each ECOD domain type |
| `pseudo-domain-map.tsv` | Pseudo-domain cluster mappings: PD cluster → target protein, residue range, TM-scores |
| `sequence_modularity.tsv` | Pairwise protein alignments with modularity category annotations |
| `structures/` | AlphaFold3-predicted structures in gzip-compressed PDB format (`.pdb.gz`) |

### RBP filtering criteria

Proteins are included if they have a non-empty `RBP-class` that is not `likely false positive`, yielding **382 proteins**.

### Domain overlap resolution

When ECOD hits overlap on the same protein:
1. Hits are grouped into non-overlapping components by coordinate overlap
2. Non-overlapping hits are kept as-is
3. For overlapping components, up to 3 hits are retained, ranked by target coverage (`tCov`) then e-value
4. Segments covered by multiple domains are labelled `multiple domains`
5. Uncovered segments are labelled `undetected`

---

## Usage

### Requirements

- Python 3 with **pandas** (`pip install pandas` or via conda)
- A modern web browser (Chrome, Firefox, Safari)

> The build script reads all TSV files and injects the processed data as JSON directly into `index.html`, so the final output is a fully self-contained single-page app — no server-side code needed at runtime.

### Running locally

```bash
git clone https://github.com/VyshakhRP/klebsiella-rbp-atlas.git
cd klebsiella-rbp-atlas

# Start a local HTTP server (required for NGL to load .pdb.gz files)
python3 -m http.server 8080
# Open http://localhost:8080
```

### Rebuilding after data changes

```bash
# With a standard Python environment that has pandas:
python3 build.py

# If using the foldseek conda environment:
conda run -n foldseek python3 build.py
```

Rebuilding regenerates `index.html` (~1.2 MB) from all TSV sources. Refresh your browser after rebuilding.

---

## Navigation quick reference

| Action | Result |
|---|---|
| Click an icon card on the home page | Opens the corresponding explorer view |
| Click the atlas title (top-left) | Returns to the home page |
| Click a protein row | Opens the detail panel with structure, sequence strip, and metadata |
| Hover over the sequence strip | Highlights the hovered residue in the 3D viewer |
| Drag the detail panel's left edge | Resizes the detail panel width |
| Drag the NGL viewport's bottom edge | Resizes the 3D viewer height |
| Click ⛶ (fullscreen) | Opens NGL in a full-screen overlay |
| Click ↙ .pdb.gz | Downloads the structure file |
| Click ↙ Export TSV | Downloads the current view as TSV |
| Click a node in the modularity network | Highlights its connections and shows protein-pair table |
| Click an edge in the modularity network | Shows protein pairs for that class-pair + category |
| Search box | Searches proteins (by ID, phage, capsule, genus, …), classes, clusters, PD clusters, ECOD domains |

---

## Libraries

| Library | Version | Purpose |
|---|---|---|
| [NGL](https://github.com/nglviewer/ngl) | 0.10.4 | 3D protein structure rendering |
| [Chart.js](https://www.chartjs.org/) | 4.4.3 | Donut charts on home page |
| [Cytoscape.js](https://cytoscape.org/) | 3.28.1 | Sequence modularity network |

All libraries are loaded from CDN; no local installation required.

---

## Citation

If you use this atlas in your research, please cite:

> Rajachandra Panicker V, et al. *A structural atlas of Klebsiella phage receptor-binding proteins reveals multi-scale modularity underlies host-range diversification.* (manuscript in preparation)

---

## Contact

Vyshakh Rajachandra Panicker — vyshakh.rajachandrapanicker@doctoral.uj.edu.pl

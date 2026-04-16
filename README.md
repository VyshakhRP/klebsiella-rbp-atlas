# Klebsiella phage RBP Atlas

An interactive browser-based atlas for exploring the structural and functional diversity of **receptor-binding proteins (RBPs)** from *Klebsiella* bacteriophages. The atlas integrates predicted protein structures, ECOD domain annotations, and phage genome metadata into a single interactive visualization.

---

## Overview

RBPs are the host-recognition modules of bacteriophages, determining which bacterial capsule serotypes a phage can infect. This atlas presents **382 RBPs** from *Klebsiella* phages, classified into **39 structural classes** and **141 clusters**, with domain architecture visualized using ECOD (Evolutionary Classification of Protein Domains) annotations and 3D structures rendered interactively in the browser.

---

## Features

- **Class-based navigation** — Browse RBPs grouped by structural class, with clusters as sub-groups within each class
- **Domain architecture visualization** — Linear domain bar diagrams colored by ECOD domain type, proportional to protein length
- **Interactive 3D structure viewer** — NGL-powered in-browser viewer with domain-matched coloring, spin, reset, and surface overlay controls
- **Host range & morphology filters** — Filter proteins by host range (Narrow / Intermediate–Broad / Unknown) and phage morphology
- **Protein & genome metadata** — Functional annotations (Pfam, ECOD, PHROGS), prediction scores, phage taxonomy, genome size, capsule type, and PubMed links
- **Search** — Find proteins by ID, genome name, or phage name; search RBP class names directly
- **Home navigation** — Click the title to return to the summary statistics page at any time

---

## Data

| File | Description |
|---|---|
| `rbp_table_v16122025.tsv` | RBP metadata: functional predictions, model quality, cluster/class assignments |
| `ecod-map.tsv` | ECOD domain hits mapped to each RBP (query = RBP, target = ECOD domain) |
| `genome_table_v25122025.tsv` | Phage genome metadata: taxonomy, morphology, host range, capsule, publications |
| `ecod-color-map.txt` | Color assignments for each ECOD domain type |
| `structures/` | AlphaFold3-predicted structures in gzip-compressed PDB format (`.pdb.gz`) |

### RBP filtering criteria
Proteins included in the atlas meet all of the following:
- Non-empty `RBP-class` annotation
- `RBP-class` is not `likely false positive`

This yields **382 proteins** from the full dataset.

### Domain overlap resolution
When multiple ECOD hits overlap on the same protein, the following logic is applied:
1. Hits are grouped into non-overlapping components based on coordinate overlap
2. Non-overlapping hits are retained as-is
3. For overlapping hits within a component, up to 3 hits are kept, prioritised by target coverage (`tCov`) and e-value
4. Protein segments covered by multiple domains are labelled `multiple domains`
5. Uncovered segments are labelled `undetected`

---

## Usage

### Requirements
- Python 3 (standard library only — no extra packages needed to serve)
- A modern web browser (Chrome, Firefox, Safari)

### Running locally

```bash
# Clone the repository
git clone https://github.com/VyshakhRP/klebsiella-phage-rbp-atlas.git
cd klebsiella-phage-rbp-atlas

# Start a local HTTP server
python3 -m http.server 8765

# Open in browser
# http://localhost:8765
```

> **Why a local server?** Browsers block direct filesystem access for security. The server makes the compressed PDB files accessible to the NGL structure viewer.

### Rebuilding the webpage

If you update any of the TSV files or the color map, regenerate `index.html` with:

```bash
pip install pandas
python3 build.py
```

Then refresh your browser.

---

## Navigating the Atlas

| Action | Result |
|---|---|
| Click an RBP class in the sidebar | Opens the class view showing all proteins grouped by cluster |
| Click the title (top-left) | Returns to the home summary page |
| Hover over a domain segment | Tooltip shows domain name and residue range |
| Click a protein row | Opens the detail panel with domain bar, 3D structure, and metadata |
| Filter buttons (top of class view) | Filter by host range or phage morphology |
| Search box | Search by protein ID, genome name, phage name, or class name |
| Spin / Reset / Surface buttons | 3D structure controls in the detail panel |
| Protein / Genome tabs | Switch between protein annotation and phage genome metadata |

---

## Domain Color Scheme

Each ECOD domain is assigned a distinct color for consistent visualization across both the 2D domain bar and the 3D structure:

| Domain | Color |
|---|---|
| Pectin lyase-like | Yellow |
| Galactose-binding domain-like | Magenta |
| Phage tail fiber protein trimerization domain | Orange |
| SGNH hydrolase | Deep blue |
| Intramolecular chaperone domain in virus tail spike | Teal |
| Ig-like domain in tailspike protein | Medium purple |
| Immunoglobulin/Fibronectin type III | Light purple |
| Concanavalin A-like lectins/glucanases | Coral pink |
| TNF-like | Dark red |
| 7-bladed beta-propeller | Beige |
| 6-bladed beta-propeller | Bisque |
| Domain in virus attachment proteins | Cyan |
| gp11/gp12 receptor-binding domain | Sky blue |
| Phage T4 gp12 N-terminal repeating units | Brown |
| Agglutinin HPA-like | Olive brown |
| gp9 C-terminal domain-related | Dark gray |
| Undetected | White |
| Multiple domains | Black |

---

## Structure Files

3D structures are AlphaFold3 predictions stored as gzip-compressed PDB files (`structures/*.pdb.gz`). The NGL viewer loads and decompresses these files on demand in the browser — no pre-processing required.

---

## Citation

If you use this atlas in your research, please cite:

> Rajachandra Panicker V, et al. *A structural atlas of Klebsiella phage receptor-binding proteins reveals multi-scale modularity underlies host-range diversification.* (manuscript submitted)

---

## Contact

Vyshakh Rajachandra Panicker (vyshakh.rajachandrapanicker@doctoral.uj.edu.pl) 


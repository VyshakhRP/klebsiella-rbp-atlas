#!/usr/bin/env python3
"""
RBP Atlas – build.py
Reads the three TSV files + color map and writes index.html.

Usage:
    python build.py
Then serve locally:
    python -m http.server 8000
    # open http://localhost:8000
"""

import ast, json, re, sys
from collections import defaultdict
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("pandas required:  pip install pandas")

# ── paths ──────────────────────────────────────────────────────────────────────
HERE         = Path(__file__).parent
RBP_TABLE    = HERE / "rbp_table_v16122025.tsv"
ECOD_MAP     = HERE / "ecod-map.tsv"
GENOME_TABLE = HERE / "genome_table_v25122025.tsv"
COLOR_FILE   = HERE / "ecod-color-map.txt"
STRUCTURES   = HERE / "structures"
OUT          = HERE / "index.html"

# ── 1. Color map ───────────────────────────────────────────────────────────────
def _rgb_hex(r, g, b):
    return "#{:02x}{:02x}{:02x}".format(int(round(r*255)), int(round(g*255)), int(round(b*255)))

raw = COLOR_FILE.read_text()
m   = re.search(r'\{(.+?)\}', raw, re.DOTALL)
if not m:
    sys.exit("Cannot parse ecod-color-map.txt")
color_map: dict = {k: _rgb_hex(*v) for k, v in ast.literal_eval("{"+m.group(1)+"}").items()}
color_map.setdefault("undetected", "#ffffff")
color_map.setdefault("multiple domains", "#000000")

# ── 2. Load tables ─────────────────────────────────────────────────────────────
rbp_df    = pd.read_csv(RBP_TABLE,    sep="\t", low_memory=False)
ecod_df   = pd.read_csv(ECOD_MAP,     sep="\t", low_memory=False)
genome_df = pd.read_csv(GENOME_TABLE, sep="\t", low_memory=False)

# ── 3. Filter to valid RBPs ────────────────────────────────────────────────────
valid_rbp = rbp_df[
    rbp_df["RBP-class"].notna() &
    (rbp_df["RBP-class"].str.strip() != "") &
    (rbp_df["RBP-class"] != "likely false positive")
].copy()
valid_ids = set(valid_rbp["proteinID"])
print(f"Valid RBPs: {len(valid_ids)}")

# ── 4. Overlap-resolution (user-supplied logic) ────────────────────────────────
df = ecod_df.copy()
df["_rank2"] = df["evalue"].astype(float)
df["qstart"] = df["qstart"].astype(int)
df["qend"]   = df["qend"].astype(int)
df["tCov"]   = df["tCov"].astype(float)

kept = []
for q, sub in df.groupby("query", sort=False):
    sub = (sub.sort_values(["qstart","qend"])
              .reset_index(drop=False)
              .rename(columns={"index":"_orig_idx"}))
    comp_ids, cur_end, comp_id = [], None, -1
    for _, r in sub.iterrows():
        if cur_end is None or r["qstart"] > cur_end:
            comp_id += 1; cur_end = r["qend"]
        else:
            cur_end = max(cur_end, r["qend"])
        comp_ids.append(comp_id)
    sub["_comp"] = comp_ids
    comp_counts = sub["_comp"].value_counts()
    sub.loc[sub["_comp"].isin(comp_counts[comp_counts == 1].index), "_status"] = "non_overlap"
    top3 = (sub[sub["_comp"].isin(comp_counts[comp_counts > 1].index)]
              .sort_values(["tCov","_rank2"], ascending=[False, True])
              .groupby("_comp", as_index=False, sort=False)
              .head(3))
    sub.loc[top3.index, "_status"] = "overlap_top3"
    kept.append(sub)

resolved = (pd.concat(kept).drop(columns=["_rank2","_comp"]).reset_index(drop=True))
resolved  = resolved[resolved["_status"].notna()].reset_index(drop=True)

# ── 5. Build per-protein tile data ─────────────────────────────────────────────
def build_tiles(sub: pd.DataFrame, qlength: int) -> list:
    borders = sorted(set(
        [1, qlength]
        + sub["qstart"].astype(int).tolist()
        + sub["qend"].astype(int).tolist()
    ))
    tiles = []
    for left, right in zip(borders[:-1], borders[1:]):
        covering = sub[(sub["qstart"] <= left) & (sub["qend"] >= right)]
        if covering.empty:
            dk, label = "undetected", "Undetected"
        else:
            dvals = covering["domain"].dropna().unique().tolist()
            if len(dvals) == 0:
                dk, label = "undetected", "Undetected"
            elif len(dvals) == 1:
                dk    = str(dvals[0])
                label = dk.split(":",1)[1] if ":" in dk else dk
            else:
                dk    = "multiple domains"
                label = " & ".join(str(d).split(":",1)[-1] for d in dvals)
        left_edge  = (left  == 1)      or (left  in set(sub["qstart"].tolist()))
        right_edge = (right == qlength) or (right in set(sub["qend"].tolist()))
        ps = left  if left_edge  else left  + 1
        pe = right if right_edge else right - 1
        if ps > pe:
            continue
        tiles.append({"s": int(ps), "e": int(pe),
                      "dk": dk,
                      "c": color_map.get(dk, color_map["undetected"]),
                      "l": label})
    return tiles

protein_tiles = {}
for q, sub in resolved[resolved["query"].isin(valid_ids)].groupby("query", sort=False):
    lengths = sub["qLen"].dropna().unique()
    qlen    = int(lengths[0]) if len(lengths) else int(sub["qend"].max())
    protein_tiles[q] = {"len": qlen, "tiles": build_tiles(sub, qlen)}

# Proteins with no ECOD hits → single undetected tile
for pid in valid_ids:
    if pid not in protein_tiles:
        row  = valid_rbp[valid_rbp["proteinID"] == pid].iloc[0]
        plen = int(row["protein_len"])
        protein_tiles[pid] = {
            "len": plen,
            "tiles": [{"s":1, "e":plen, "dk":"undetected",
                       "c": color_map["undetected"], "l":"Undetected"}]
        }

# ── 6. Genome metadata ─────────────────────────────────────────────────────────
valid_rbp = valid_rbp.copy()
valid_rbp["Genome_name"] = valid_rbp["proteinID"].str.replace(r"(_PROTEIN_.*|_gp.*)", "", regex=True)
genome_lookup = genome_df.set_index("Genome_name").to_dict("index")

def hr_label(raw):
    if str(raw) in ("Intermediate","Broad"): return "I/B"
    if str(raw) == "Narrow": return "Narrow"
    return "Unknown"

# ── 7. Build class → cluster → protein hierarchy ───────────────────────────────
PMETA = ["protein_len","RBP","PFAM_function",
         "ECOD_function","PHROGS1_function","PHROGS2_function","preds","score",
         "ranking_score","RBPseg(YES/NO)","Avg_MergeScore","Model_Quality","TF-class"]

hierarchy = defaultdict(lambda: defaultdict(list))

for _, row in valid_rbp.iterrows():
    pid     = row["proteinID"]
    cls     = row["RBP-class"]
    cluster = str(row.get("RBP-cluster","") or "Unknown")
    gname   = row["Genome_name"]
    gm      = genome_lookup.get(gname, {})
    hr_raw  = str(gm.get("Host_range",""))
    hr      = hr_label(hr_raw)
    has_pdb = (STRUCTURES / f"{pid}.pdb.gz").exists()

    pmeta = {}
    for c in PMETA:
        if c in row.index:
            v = row[c]
            pmeta[c.replace("(YES/NO)","")] = "" if pd.isna(v) else str(v)

    gmeta = {
        "Phage":         str(gm.get("Phage","")),
        "Genome_name":   str(gm.get("Genome_name", gname)),
        "Genome_ID":     str(gm.get("Genome_ID","")),
        "Genus":         str(gm.get("Genus","")),
        "Phage_cluster": str(gm.get("Phage_cluster","")),
        "Genome_size_bp":str(gm.get("Genome_size_(bp)","")),
        "Morphology":    str(gm.get("Morphology","")),
        "Capsule":       str(gm.get("Capsule","")),
        "Observation":   str(gm.get("Observation","")),
        "PMID":          str(gm.get("Publication(PMID)","")),
        "Host_range":    hr_raw,
        "Halo":          str(gm.get("Halo","")),
        "HR":            hr,
    }

    ti = protein_tiles[pid]
    hierarchy[cls][cluster].append({
        "id":   pid,
        "gn":    gname,
        "plen":  int(row["protein_len"]),
        "pdb":   f"structures/{pid}.pdb.gz" if has_pdb else None,
        "hr":    hr,
        "morph": gmeta["Morphology"],
        "clid":  cluster,
        "tiles": ti["tiles"],
        "tlen":  ti["len"],
        "pm":    pmeta,
        "gm":    gmeta,
    })

# ── 8. Serialise ───────────────────────────────────────────────────────────────
def _ckey(k):
    m2 = re.match(r"RBP(\d+)$", k)
    return int(m2.group(1)) if m2 else 9999

classes_list = []
for cls in sorted(hierarchy):
    clusters = hierarchy[cls]
    cl = [{"id":cid, "proteins":sorted(clusters[cid], key=lambda p: p["id"])}
          for cid in sorted(clusters, key=_ckey)]
    n  = sum(len(c["proteins"]) for c in cl)
    classes_list.append({"name": cls, "clusters": cl, "n": n})

atlas_data = {"colorMap": color_map, "classes": classes_list}
data_json  = json.dumps(atlas_data, separators=(",",":"), ensure_ascii=False)
print(f"JSON: {len(data_json)//1024} KB  |  {len(valid_ids)} proteins  |  {len(classes_list)} classes")

# ── 9. HTML template ───────────────────────────────────────────────────────────
HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Klebsiella phage RBP Atlas</title>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{display:flex;height:100vh;overflow:hidden;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;font-size:13px;background:#f4f6f9;color:#222}

/* ── Sidebar ── */
#sidebar{width:270px;min-width:270px;background:#1a2b3c;color:#c8d8e8;display:flex;flex-direction:column;overflow:hidden;z-index:10}
.sb-header{padding:18px 16px 12px;border-bottom:1px solid #2d4157}
.sb-header h1{font-size:18px;font-weight:700;color:#fff;letter-spacing:.3px}
.sb-header .sub{font-size:11px;color:#7a9ab8;margin-top:3px}
.sb-search{padding:10px 12px;border-bottom:1px solid #2d4157}
.sb-search input{width:100%;padding:7px 10px;border-radius:5px;border:none;background:#253d52;color:#e0eaf4;font-size:12px;outline:none}
.sb-search input::placeholder{color:#5a7a94}
.sb-search input:focus{background:#2d4d66;box-shadow:0 0 0 2px #4a90d9}
#class-nav{flex:1;overflow-y:auto;padding:6px 0}
#class-nav::-webkit-scrollbar{width:4px}
#class-nav::-webkit-scrollbar-thumb{background:#2d4157;border-radius:2px}
.cls-item{cursor:pointer}
.cls-label{display:flex;align-items:center;padding:8px 14px;gap:8px;transition:background .15s}
.cls-label:hover{background:#253d52}
.cls-label.active{background:#2d5a8e;color:#fff}
.cls-arrow{font-size:10px;color:#5a7a94;transition:transform .2s;flex-shrink:0}
.cls-arrow.open{transform:rotate(90deg)}
.cls-name{flex:1;font-size:12px;line-height:1.4;color:#c8d8e8}
.cls-label.active .cls-name{color:#fff}
.cls-count{background:#2d4157;color:#7a9ab8;font-size:10px;font-weight:600;padding:1px 6px;border-radius:10px;flex-shrink:0}
.cls-label.active .cls-count{background:#1e4070;color:#9ac0f0}
.cluster-list{display:none;background:#152333}
.cluster-list.open{display:block}
.cl-item{display:flex;align-items:center;padding:5px 14px 5px 30px;cursor:pointer;gap:6px;transition:background .12s}
.cl-item:hover{background:#1e3347}
.cl-item.active{background:#1e4070}
.cl-dot{width:7px;height:7px;border-radius:50%;background:#4a90d9;flex-shrink:0}
.cl-id{font-size:11px;font-weight:600;color:#9ac0f0}
.cl-n{font-size:10px;color:#5a7a94;margin-left:auto}

/* ── Main panel ── */
#main{flex:1;overflow-y:auto;padding:0;display:flex;flex-direction:column;min-width:0}
#main::-webkit-scrollbar{width:6px}
#main::-webkit-scrollbar-thumb{background:#c8d0da;border-radius:3px}

/* Welcome */
#welcome{display:flex;flex-direction:column;align-items:center;justify-content:center;height:100%;color:#666;gap:12px}
#welcome .big{font-size:28px;font-weight:700;color:#1a2b3c}
#welcome .hint{font-size:14px;color:#888;max-width:380px;text-align:center;line-height:1.6}
.welcome-stats{display:flex;gap:24px;margin-top:8px}
.ws-card{background:#fff;border:1px solid #e0e6ec;border-radius:8px;padding:14px 22px;text-align:center}
.ws-num{font-size:28px;font-weight:700;color:#2d5a8e}
.ws-lbl{font-size:11px;color:#888;margin-top:2px}

/* Class view */
#class-view{display:none;flex-direction:column;height:100%}
.cv-header{background:#fff;border-bottom:1px solid #e0e6ec;padding:14px 20px;display:flex;align-items:center;gap:16px;flex-wrap:wrap;flex-shrink:0}
.cv-title{font-size:16px;font-weight:700;color:#1a2b3c;flex:1}
.cv-meta{font-size:12px;color:#888}
.filter-row{display:flex;gap:8px;align-items:center;flex-wrap:wrap}
.filter-label{font-size:11px;color:#888;font-weight:600}
.filter-btn{font-size:11px;padding:3px 10px;border-radius:12px;border:1px solid #dde3ea;background:#f8f9fa;cursor:pointer;transition:all .15s;color:#555}
.filter-btn:hover{border-color:#4a90d9;color:#2d5a8e}
.filter-btn.active{background:#2d5a8e;color:#fff;border-color:#2d5a8e}
.cv-body{flex:1;overflow-y:auto;padding:16px 20px}
.cv-body::-webkit-scrollbar{width:6px}
.cv-body::-webkit-scrollbar-thumb{background:#c8d0da;border-radius:3px}

/* Cluster section */
.cluster-section{margin-bottom:18px}
.cluster-header{display:flex;align-items:center;gap:10px;padding:8px 12px;background:#fff;border:1px solid #e0e6ec;border-radius:6px 6px 0 0;border-bottom:2px solid #2d5a8e}
.cluster-id{font-size:13px;font-weight:700;color:#2d5a8e}
.cluster-n{font-size:11px;color:#888}
.cluster-proteins{background:#fff;border:1px solid #e0e6ec;border-top:none;border-radius:0 0 6px 6px;overflow:hidden}
.prot-row{display:flex;align-items:center;gap:10px;padding:6px 12px;border-bottom:1px solid #f0f3f7;cursor:pointer;transition:background .12s}
.prot-row:last-child{border-bottom:none}
.prot-row:hover{background:#f0f4ff}
.prot-row.selected{background:#e8f0fe;border-left:3px solid #4a90d9}
.prot-id{font-family:'SF Mono',Monaco,monospace;font-size:11px;color:#333;min-width:210px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.bar-col{flex:1;min-width:0}
.hr-badge{font-size:10px;font-weight:700;padding:2px 7px;border-radius:10px;flex-shrink:0;min-width:46px;text-align:center}
.hr-Narrow{background:#fff0f0;color:#c0392b}
.hr-IB{background:#f0fff4;color:#27ae60}
.hr-Unknown{background:#f5f5f5;color:#999}
.morph-tag{font-size:10px;color:#888;min-width:72px;text-align:right;flex-shrink:0}
.no-pdb{font-size:10px;color:#ccc}

/* Domain bar */
.domain-bar{display:flex;height:16px;border-radius:3px;overflow:hidden;border:1px solid #dde3ea;min-width:0}
.domain-bar.large{height:26px;border-radius:4px}
.dseg{height:100%;cursor:default;transition:filter .1s}
.dseg:hover{filter:brightness(0.88)}
.dseg.undetected{border-right:1px solid #ececec}

/* Detail panel */
#detail{width:440px;min-width:440px;background:#fff;border-left:1px solid #dde3ea;display:flex;flex-direction:column;transform:translateX(100%);transition:transform .28s ease;position:fixed;right:0;top:0;bottom:0;z-index:200;overflow:hidden}
#detail.open{transform:translateX(0)}
.det-header{padding:14px 16px 12px;border-bottom:1px solid #eef1f6;flex-shrink:0;background:#fff}
.det-close{float:right;background:none;border:none;cursor:pointer;font-size:18px;color:#aaa;line-height:1;padding:2px 6px;border-radius:4px}
.det-close:hover{background:#f0f3f7;color:#555}
.det-pid{font-family:'SF Mono',Monaco,monospace;font-size:13px;font-weight:700;color:#1a2b3c;margin-bottom:4px;padding-right:30px}
.det-tags{display:flex;gap:6px;flex-wrap:wrap}
.det-tag{font-size:10px;padding:2px 8px;border-radius:10px;background:#f0f3f7;color:#555}
.det-tag.cls{background:#e8f0fe;color:#2d5a8e;font-weight:600}
.det-body{flex:1;overflow-y:auto;padding:0}
.det-body::-webkit-scrollbar{width:5px}
.det-body::-webkit-scrollbar-thumb{background:#dde3ea;border-radius:3px}
.det-section{padding:14px 16px;border-bottom:1px solid #eef1f6}
.det-section h3{font-size:12px;font-weight:700;color:#888;text-transform:uppercase;letter-spacing:.6px;margin-bottom:10px}
.bar-ruler{display:flex;justify-content:space-between;font-size:10px;color:#aaa;margin-top:3px;font-family:monospace}
.domain-legend{display:flex;flex-wrap:wrap;gap:8px;margin-top:10px}
.leg-item{display:flex;align-items:center;gap:5px;font-size:11px;color:#444}
.leg-swatch{width:14px;height:14px;border-radius:2px;border:1px solid rgba(0,0,0,.12);flex-shrink:0}

/* NGL */
#ngl-viewport{height:280px;background:#f8f9fa;border-radius:6px;overflow:hidden;border:1px solid #e0e6ec;position:relative}
.ngl-loading{position:absolute;inset:0;display:flex;align-items:center;justify-content:center;font-size:12px;color:#aaa}
.ngl-controls{display:flex;gap:8px;margin-top:8px;align-items:center}
.ngl-btn{font-size:11px;padding:4px 12px;border-radius:4px;border:1px solid #dde3ea;background:#f8f9fa;cursor:pointer;transition:all .15s;color:#555}
.ngl-btn:hover{border-color:#4a90d9;color:#2d5a8e}
.ngl-btn.active{background:#2d5a8e;color:#fff;border-color:#2d5a8e}
.ngl-no-struct{display:flex;align-items:center;justify-content:center;height:100%;color:#bbb;font-size:12px}

/* Metadata tabs */
.tabs{display:flex;border-bottom:2px solid #eef1f6;margin-bottom:12px}
.tab-btn{padding:7px 16px;font-size:12px;font-weight:600;color:#888;border:none;background:none;cursor:pointer;border-bottom:2px solid transparent;margin-bottom:-2px;transition:all .15s}
.tab-btn:hover{color:#2d5a8e}
.tab-btn.active{color:#2d5a8e;border-bottom-color:#2d5a8e}
.tab-content{display:none}
.tab-content.active{display:block}
.meta-table{width:100%;border-collapse:collapse;font-size:12px}
.meta-table tr:not(:last-child) td{border-bottom:1px solid #f0f3f7}
.meta-table td{padding:5px 8px;vertical-align:top}
.meta-table td:first-child{color:#888;font-weight:600;white-space:nowrap;width:140px}
.meta-table td:last-child{color:#333;word-break:break-word}
.pmid-link{color:#4a90d9;text-decoration:none}
.pmid-link:hover{text-decoration:underline}

/* Tooltip */
#tooltip{position:fixed;pointer-events:none;background:rgba(20,30,45,.9);color:#e8f0f8;font-size:11px;padding:6px 10px;border-radius:5px;max-width:260px;line-height:1.5;z-index:9999;display:none;box-shadow:0 2px 8px rgba(0,0,0,.3)}

/* Search results */
#search-results{position:absolute;top:100%;left:0;right:0;background:#1a2b3c;border:1px solid #2d4157;border-radius:0 0 6px 6px;max-height:280px;overflow-y:auto;z-index:100;display:none}
.sb-search{position:relative}
.sr-item{padding:7px 14px;cursor:pointer;border-bottom:1px solid #1e3347;transition:background .12s}
.sr-item:hover{background:#253d52}
.sr-pid{font-family:monospace;font-size:11px;color:#9ac0f0}
.sr-cls{font-size:10px;color:#5a7a94;margin-top:1px}
</style>
</head>
<body>

<aside id="sidebar">
  <div class="sb-header">
    <h1 id="site-title" style="cursor:pointer" onclick="goHome()">Klebsiella phage RBP Atlas</h1>
    <div class="sub" id="sb-stats">Loading…</div>
  </div>
  <div class="sb-search">
    <input id="search-input" type="text" placeholder="Search protein or class…" autocomplete="off">
    <div id="search-results"></div>
  </div>
  <nav id="class-nav"></nav>
</aside>

<main id="main">
  <div id="welcome">
    <div class="big">Klebsiella phage RBP Atlas</div>
    <div class="hint">Select an RBP class from the left panel to explore domain architectures and 3D structures.</div>
    <div class="welcome-stats" id="w-stats"></div>
  </div>
  <div id="class-view"></div>
</main>

<aside id="detail">
  <div class="det-header">
    <button class="det-close" id="close-detail">✕</button>
    <div class="det-pid" id="det-pid"></div>
    <div class="det-tags" id="det-tags"></div>
  </div>
  <div class="det-body">
    <div class="det-section" id="det-arch-sec">
      <h3>Domain Architecture</h3>
      <div id="det-bar"></div>
      <div id="det-ruler" class="bar-ruler"></div>
      <div id="det-legend" class="domain-legend"></div>
    </div>
    <div class="det-section" id="det-struct-sec">
      <h3>3D Structure</h3>
      <div id="ngl-viewport"><div class="ngl-loading">Loading structure…</div></div>
      <div class="ngl-controls">
        <button class="ngl-btn" id="btn-spin">⟳ Spin</button>
        <button class="ngl-btn" id="btn-reset">⊕ Reset</button>
        <button class="ngl-btn" id="btn-surf">◼ Surface</button>
      </div>
    </div>
    <div class="det-section">
      <div class="tabs">
        <button class="tab-btn active" data-tab="protein">Protein</button>
        <button class="tab-btn" data-tab="genome">Genome</button>
      </div>
      <div id="tab-protein" class="tab-content active"></div>
      <div id="tab-genome"  class="tab-content"></div>
    </div>
  </div>
</aside>

<div id="tooltip"></div>

<script src="https://cdn.jsdelivr.net/npm/ngl@0.10.4/dist/ngl.min.js"></script>
<script>
const D = __ATLAS_DATA__;

// ── State ──────────────────────────────────────────────────────────────────────
const S = {
  cls: null,       // selected class name
  prot: null,      // selected protein id
  filters: { hr: "all", morph: "all" },
  nglStage: null,
  nglComp: null,
  spinning: false,
  surfaceRep: null,
};

// ── Lookup helpers ────────────────────────────────────────────────────────────
const classMap   = new Map(D.classes.map(c => [c.name, c]));
const proteinMap = new Map();
D.classes.forEach(c => c.clusters.forEach(cl => cl.proteins.forEach(p => proteinMap.set(p.id, p))));

// ── Init ──────────────────────────────────────────────────────────────────────
document.addEventListener("DOMContentLoaded", () => {
  const total = D.classes.reduce((a, c) => a + c.n, 0);
  document.getElementById("sb-stats").textContent = `${total} proteins · ${D.classes.length} classes`;

  // welcome stats
  const totalClusters = D.classes.reduce((a, c) => a + c.clusters.length, 0);
  const ecodDomains = new Set();
  D.classes.forEach(c => c.clusters.forEach(cl => cl.proteins.forEach(p =>
    p.tiles.forEach(t => { if (t.dk !== "undetected" && t.dk !== "multiple domains") ecodDomains.add(t.dk); })
  )));
  const wStats = document.getElementById("w-stats");
  [["RBPs", total], ["RBP classes", D.classes.length], ["RBP clusters", totalClusters], ["ECOD domains", ecodDomains.size]]
    .forEach(([l,n]) => {
      wStats.innerHTML += `<div class="ws-card"><div class="ws-num">${n}</div><div class="ws-lbl">${l}</div></div>`;
    });

  renderSidebar();
  setupSearch();
  setupDetailButtons();
  setupTooltip();
});

// ── Sidebar ───────────────────────────────────────────────────────────────────
function renderSidebar() {
  const nav = document.getElementById("class-nav");
  nav.innerHTML = D.classes.map(c => `
    <div class="cls-item" data-cls="${esc(c.name)}">
      <div class="cls-label" onclick="selectClass('${esc(c.name)}')">
        <span class="cls-arrow">▶</span>
        <span class="cls-name">${esc(c.name)}</span>
        <span class="cls-count">${c.n}</span>
      </div>
      <div class="cluster-list">
        ${c.clusters.map(cl => `
          <div class="cl-item" onclick="scrollToCluster('${esc(cl.id)}')" data-clid="${esc(cl.id)}">
            <span class="cl-dot"></span>
            <span class="cl-id">${esc(cl.id)}</span>
            <span class="cl-n">${cl.proteins.length}</span>
          </div>`).join("")}
      </div>
    </div>`).join("");
}

function selectClass(name) {
  if (S.cls === name) return;
  S.cls = name;
  S.filters = { hr: "all", morph: "all" };

  // sidebar active + expand
  document.querySelectorAll(".cls-label").forEach(el => {
    const active = el.closest(".cls-item").dataset.cls === name;
    el.classList.toggle("active", active);
    el.querySelector(".cls-arrow").classList.toggle("open", active);
    el.nextElementSibling.classList.toggle("open", active);
  });

  renderClassView(name);
}

function scrollToCluster(cid) {
  const el = document.getElementById("cl-" + cid);
  if (el) el.scrollIntoView({behavior:"smooth", block:"start"});
  // highlight active cluster link
  document.querySelectorAll(".cl-item").forEach(e =>
    e.classList.toggle("active", e.dataset.clid === cid));
}

// ── Class view ────────────────────────────────────────────────────────────────
function renderClassView(name) {
  const cv = document.getElementById("class-view");
  document.getElementById("welcome").style.display = "none";
  cv.style.display = "flex";
  cv.style.flexDirection = "column";
  cv.style.height = "100%";

  const cls = classMap.get(name);
  const morphs = [...new Set(cls.clusters.flatMap(cl =>
    cl.proteins.map(p => p.morph).filter(Boolean)))].sort();

  cv.innerHTML = `
    <div class="cv-header">
      <div>
        <div class="cv-title">${esc(name)}</div>
        <div class="cv-meta">${cls.clusters.length} cluster${cls.clusters.length>1?"s":""} · ${cls.n} protein${cls.n>1?"s":""}</div>
      </div>
      <div class="filter-row">
        <span class="filter-label">Host range:</span>
        ${["all","Narrow","I/B","Unknown"].map(v =>
          `<button class="filter-btn${S.filters.hr===v?" active":""}" onclick="setFilter('hr','${v}')">${v==="all"?"All":v}</button>`
        ).join("")}
        <span class="filter-label" style="margin-left:8px">Morphology:</span>
        <button class="filter-btn${S.filters.morph==="all"?" active":""}" onclick="setFilter('morph','all')">All</button>
        ${morphs.map(m =>
          `<button class="filter-btn${S.filters.morph===m?" active":""}" onclick="setFilter('morph','${esc(m)}')">${esc(m)}</button>`
        ).join("")}
      </div>
    </div>
    <div class="cv-body" id="cv-body"></div>`;

  renderClusters(cls);
}

function renderClusters(cls) {
  const body = document.getElementById("cv-body");
  if (!body) return;
  const { hr, morph } = S.filters;

  let html = "";
  let visible = 0;
  cls.clusters.forEach(cl => {
    const prots = cl.proteins.filter(p =>
      (hr === "all" || p.hr === hr) &&
      (morph === "all" || p.morph === morph));
    if (!prots.length) return;
    visible += prots.length;
    html += `
      <div class="cluster-section" id="cl-${esc(cl.id)}">
        <div class="cluster-header">
          <span class="cluster-id">${esc(cl.id)}</span>
          <span class="cluster-n">${prots.length} protein${prots.length>1?"s":""}</span>
        </div>
        <div class="cluster-proteins">
          ${prots.map(p => proteinRowHTML(p)).join("")}
        </div>
      </div>`;
  });

  if (!visible) {
    html = `<div style="padding:40px;text-align:center;color:#aaa">No proteins match the current filters.</div>`;
  }
  body.innerHTML = html;
}

function proteinRowHTML(p) {
  const sel = S.prot === p.id ? " selected" : "";
  const hrCls = "hr-" + p.hr.replace("/","");
  return `<div class="prot-row${sel}" onclick="selectProtein('${esc(p.id)}')" id="pr-${esc(p.id)}">
    <span class="prot-id" title="${esc(p.id)}">${esc(p.id)}</span>
    <span class="bar-col">${domainBarHTML(p, 16)}</span>
    <span class="hr-badge ${hrCls}">${esc(p.hr)}</span>
    <span class="morph-tag">${esc(p.morph||"")}</span>
  </div>`;
}

function setFilter(key, val) {
  S.filters[key] = val;
  const cls = classMap.get(S.cls);
  if (!cls) return;
  // re-render header and clusters
  renderClassView(S.cls);
}

// ── Domain bar ────────────────────────────────────────────────────────────────
function domainBarHTML(p, height) {
  const len = p.tlen || p.plen || 1;
  const segs = p.tiles.map(t => {
    const w = ((t.e - t.s + 1) / len * 100).toFixed(3);
    const border = t.dk === "undetected" ? "border-right:1px solid #ececec;" : "";
    const label  = `${t.l} (${t.s}–${t.e})`;
    return `<div class="dseg${t.dk==="undetected"?" undetected":""}"
      style="flex:${t.e-t.s+1};background:${t.c};${border}height:${height}px"
      data-tip="${esc(label)}"></div>`;
  }).join("");
  return `<div class="domain-bar" style="height:${height}px">${segs}</div>`;
}

// ── Protein detail ────────────────────────────────────────────────────────────
function selectProtein(id) {
  S.prot = id;
  const p = proteinMap.get(id);
  if (!p) return;

  // highlight row
  document.querySelectorAll(".prot-row").forEach(el =>
    el.classList.toggle("selected", el.id === "pr-"+id));

  openDetail(p);
}

function openDetail(p) {
  document.getElementById("detail").classList.add("open");

  // header
  document.getElementById("det-pid").textContent = p.id;
  const tags = document.getElementById("det-tags");
  tags.innerHTML = "";
  const addTag = (txt, cls="") =>
    txt && (tags.innerHTML += `<span class="det-tag ${cls}">${esc(String(txt))}</span>`);
  addTag(S.cls, "cls");
  addTag(p.clid);
  addTag(p.hr);
  addTag(p.morph);
  addTag(p.plen + " aa");

  // domain bar (large)
  document.getElementById("det-bar").innerHTML = domainBarHTML(p, 26);
  document.getElementById("det-ruler").innerHTML =
    `<span>1</span><span>${p.tlen}</span>`;

  // legend
  const seen = new Map();
  p.tiles.forEach(t => { if (!seen.has(t.dk)) seen.set(t.dk, t); });
  document.getElementById("det-legend").innerHTML =
    [...seen.values()].map(t =>
      `<div class="leg-item">
        <div class="leg-swatch" style="background:${t.c}"></div>
        <span>${esc(t.l)}</span>
       </div>`).join("");

  // metadata tables
  const PM_LABELS = {
    protein_len:"Length (aa)", RBP:"RBP prediction",
    PFAM_function:"Pfam function", ECOD_function:"ECOD function",
    PHROGS1_function:"PHROGS1 function", PHROGS2_function:"PHROGS2 function",
    preds:"Prediction confidence", score:"Score", ranking_score:"Ranking score",
    RBPseg:"RBPseg", Avg_MergeScore:"Avg merge score",
    Model_Quality:"Model quality", "TF-class":"TF class"
  };
  const GM_LABELS = {
    Phage:"Phage name", Genome_name:"Genome", Genome_ID:"Genome ID",
    Genus:"Genus", Phage_cluster:"Phage cluster", Genome_size_bp:"Genome size (bp)",
    Morphology:"Morphology", Capsule:"Capsule", Observation:"Observation",
    PMID:"PubMed ID", Host_range:"Host range", Halo:"Halo", HR:"HR group"
  };

  const makeTable = (obj, labels) => {
    let rows = "";
    for (const [k, lbl] of Object.entries(labels)) {
      const v = obj[k] || "";
      if (!v || v === "nan" || v === "NaN") continue;
      let vhtml = esc(v);
      if (k === "PMID" && v && v !== "nan") {
        const ids = v.split(/[,;]\s*/);
        vhtml = ids.map(id => {
          const t = id.trim();
          if (!t) return "";
          // Only linkify numeric PMIDs
          return /^\d+$/.test(t)
            ? `<a class="pmid-link" href="https://pubmed.ncbi.nlm.nih.gov/${t}" target="_blank">${t}</a>`
            : esc(t);
        }).join(", ");
      }
      rows += `<tr><td>${esc(lbl)}</td><td>${vhtml}</td></tr>`;
    }
    return `<table class="meta-table"><tbody>${rows}</tbody></table>`;
  };

  document.getElementById("tab-protein").innerHTML = makeTable(p.pm, PM_LABELS);
  document.getElementById("tab-genome").innerHTML  = makeTable(p.gm, GM_LABELS);

  // 3D structure
  initNGL(p);
}

function closeDetail() {
  document.getElementById("detail").classList.remove("open");
  S.prot = null;
  if (S.nglStage) { S.nglStage.dispose(); S.nglStage = null; }
  S.spinning = false;
  document.querySelectorAll(".prot-row").forEach(el => el.classList.remove("selected"));
}

// ── NGL viewer ────────────────────────────────────────────────────────────────
function initNGL(p) {
  const el = document.getElementById("ngl-viewport");
  if (S.nglStage) { S.nglStage.dispose(); S.nglStage = null; }
  S.spinning = false;
  document.getElementById("btn-spin").classList.remove("active");
  S.surfaceRep = null;
  document.getElementById("btn-surf").classList.remove("active");

  if (!p.pdb) {
    el.innerHTML = '<div class="ngl-no-struct">Structure not available</div>';
    return;
  }
  el.innerHTML = '<div class="ngl-loading">Loading structure…</div>';

  try {
    S.nglStage = new NGL.Stage(el, { backgroundColor: "white" });
    S.nglStage.loadFile(p.pdb, { ext: "pdb.gz" }).then(comp => {
      // Remove the loading overlay now that NGL has taken over the element
      const loading = el.querySelector(".ngl-loading");
      if (loading) loading.remove();
      S.nglComp = comp;
      const scheme = buildColorScheme(p.tiles);
      comp.addRepresentation("cartoon", { color: scheme, smoothSheet: true });
      comp.autoView();
    }).catch(err => {
      el.innerHTML = `<div class="ngl-no-struct">Error loading structure</div>`;
      console.error(err);
    });
  } catch(e) {
    el.innerHTML = `<div class="ngl-no-struct">NGL not available – serve via HTTP</div>`;
  }
}

function buildColorScheme(tiles) {
  const sel = tiles.map(t => [t.c, `${t.s}-${t.e}`]);
  sel.push(["#aaaaaa", "*"]);
  return NGL.ColormakerRegistry.addSelectionScheme(sel, "dc_" + Date.now());
}

function setupDetailButtons() {
  document.getElementById("close-detail").onclick = closeDetail;

  document.getElementById("btn-spin").onclick = function() {
    S.spinning = !S.spinning;
    this.classList.toggle("active", S.spinning);
    if (S.nglStage) S.nglStage.setSpin(S.spinning ? [0,1,0] : null);
  };

  document.getElementById("btn-reset").onclick = () => {
    if (S.nglComp) S.nglComp.autoView();
  };

  document.getElementById("btn-surf").onclick = function() {
    if (!S.nglComp) return;
    if (S.surfaceRep) {
      S.nglComp.removeRepresentation(S.surfaceRep);
      S.surfaceRep = null;
      this.classList.remove("active");
    } else {
      S.surfaceRep = S.nglComp.addRepresentation("surface",
        { opacity: 0.3, color: "white", side: "front" });
      this.classList.add("active");
    }
  };

  // tabs
  document.querySelectorAll(".tab-btn").forEach(btn => {
    btn.onclick = function() {
      document.querySelectorAll(".tab-btn").forEach(b => b.classList.remove("active"));
      document.querySelectorAll(".tab-content").forEach(c => c.classList.remove("active"));
      this.classList.add("active");
      document.getElementById("tab-"+this.dataset.tab).classList.add("active");
    };
  });
}

// ── Tooltip ───────────────────────────────────────────────────────────────────
function setupTooltip() {
  const tip = document.getElementById("tooltip");
  document.addEventListener("mouseover", e => {
    const el = e.target.closest("[data-tip]");
    if (el) {
      tip.textContent = el.dataset.tip;
      tip.style.display = "block";
    }
  });
  document.addEventListener("mousemove", e => {
    tip.style.left = (e.clientX + 14) + "px";
    tip.style.top  = (e.clientY + 14) + "px";
  });
  document.addEventListener("mouseout", e => {
    if (e.target.closest("[data-tip]")) tip.style.display = "none";
  });
}

// ── Search ────────────────────────────────────────────────────────────────────
function setupSearch() {
  const input  = document.getElementById("search-input");
  const resDiv = document.getElementById("search-results");

  input.addEventListener("input", () => {
    const q = input.value.trim().toLowerCase();
    if (q.length < 2) { resDiv.style.display = "none"; return; }

    const hits = [];
    for (const [id, p] of proteinMap) {
      if (hits.length >= 40) break;
      if (id.toLowerCase().includes(q) ||
          p.gn.toLowerCase().includes(q) ||
          (p.gm.Phage||"").toLowerCase().includes(q)) {
        hits.push(p);
      }
    }
    // also class name matches
    const clsHits = D.classes.filter(c => c.name.toLowerCase().includes(q));

    if (!hits.length && !clsHits.length) { resDiv.style.display = "none"; return; }

    let html = "";
    clsHits.slice(0,5).forEach(c => {
      html += `<div class="sr-item" onclick="selectClass('${esc(c.name)}');clearSearch()">
        <div class="sr-pid">📂 ${esc(c.name)}</div>
        <div class="sr-cls">${c.n} proteins</div></div>`;
    });
    hits.forEach(p => {
      html += `<div class="sr-item" onclick="jumpToProtein('${esc(p.id)}');clearSearch()">
        <div class="sr-pid">${esc(p.id)}</div>
        <div class="sr-cls">${esc(S.cls || classOfProtein(p.id))}</div></div>`;
    });
    resDiv.innerHTML = html;
    resDiv.style.display = "block";
  });

  document.addEventListener("click", e => {
    if (!e.target.closest(".sb-search")) resDiv.style.display = "none";
  });
}

function clearSearch() {
  document.getElementById("search-input").value = "";
  document.getElementById("search-results").style.display = "none";
}

function classOfProtein(id) {
  for (const c of D.classes) {
    for (const cl of c.clusters) {
      if (cl.proteins.some(p => p.id === id)) return c.name;
    }
  }
  return "";
}

function jumpToProtein(id) {
  const cls = classOfProtein(id);
  if (!cls) return;
  selectClass(cls);
  setTimeout(() => {
    selectProtein(id);
    const row = document.getElementById("pr-"+id);
    if (row) row.scrollIntoView({behavior:"smooth", block:"center"});
  }, 100);
}

// ── Utility ───────────────────────────────────────────────────────────────────
function goHome() {
  S.cls = null; S.prot = null;
  closeDetail();
  document.getElementById("class-view").style.display = "none";
  document.getElementById("welcome").style.display = "";
  document.querySelectorAll(".cls-label").forEach(el => {
    el.classList.remove("active");
    el.querySelector(".cls-arrow").classList.remove("open");
    el.nextElementSibling.classList.remove("open");
  });
}

function esc(s) {
  return String(s||"")
    .replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;")
    .replace(/"/g,"&quot;").replace(/'/g,"&#39;");
}
</script>
</body>
</html>
"""

# ── 10. Inject data + write ────────────────────────────────────────────────────
html_out = HTML.replace("__ATLAS_DATA__", data_json)
OUT.write_text(html_out, encoding="utf-8")
print(f"Written: {OUT}  ({OUT.stat().st_size//1024} KB)")

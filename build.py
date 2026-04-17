#!/usr/bin/env python3
"""
RBP Atlas – build.py
Usage:  conda run -n foldseek python3 build.py
Serve:  python3 -m http.server 8080  → http://localhost:8080
"""

import ast, json, re, sys
from collections import defaultdict
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("pandas required:  pip install pandas")

HERE         = Path(__file__).parent
RBP_TABLE    = HERE / "rbp_table_v16122025.tsv"
ECOD_MAP     = HERE / "ecod-map.tsv"
GENOME_TABLE = HERE / "genome_table_v25122025.tsv"
COLOR_FILE   = HERE / "ecod-color-map.txt"
PD_MAP       = HERE / "pseudo-domain-map.tsv"
SM_MAP       = HERE / "sequence_modularity.tsv"
STRUCTURES   = HERE / "structures"
OUT          = HERE / "index.html"

# ── 1. Color map ───────────────────────────────────────────────────────────────
def _rgb_hex(r, g, b):
    return "#{:02x}{:02x}{:02x}".format(int(round(r*255)), int(round(g*255)), int(round(b*255)))

raw = COLOR_FILE.read_text()
m   = re.search(r'\{(.+?)\}', raw, re.DOTALL)
if not m:
    sys.exit("Cannot parse ecod-color-map.txt")
color_map = {k: _rgb_hex(*v) for k, v in ast.literal_eval("{"+m.group(1)+"}").items()}
color_map.setdefault("undetected", "#ffffff")
color_map.setdefault("multiple domains", "#000000")

# ── 2. Load tables ─────────────────────────────────────────────────────────────
rbp_df    = pd.read_csv(RBP_TABLE,    sep="\t", low_memory=False)
ecod_df   = pd.read_csv(ECOD_MAP,     sep="\t", low_memory=False)
genome_df = pd.read_csv(GENOME_TABLE, sep="\t", low_memory=False)
pd_df     = pd.read_csv(PD_MAP,       sep="\t", low_memory=False)

# ── 3. Filter valid RBPs ───────────────────────────────────────────────────────
valid_rbp = rbp_df[
    rbp_df["RBP-class"].notna() &
    (rbp_df["RBP-class"].str.strip() != "") &
    (rbp_df["RBP-class"] != "likely false positive")
].copy()
valid_ids = set(valid_rbp["proteinID"])
print(f"Valid RBPs: {len(valid_ids)}")

# ── 4. ECOD overlap resolution ─────────────────────────────────────────────────
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
              .groupby("_comp", as_index=False, sort=False).head(3))
    sub.loc[top3.index, "_status"] = "overlap_top3"
    kept.append(sub)

resolved = pd.concat(kept).drop(columns=["_rank2","_comp"]).reset_index(drop=True)
resolved  = resolved[resolved["_status"].notna()].reset_index(drop=True)

# ── 5. Tile data ───────────────────────────────────────────────────────────────
def build_tiles(sub, qlength):
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
            if not dvals:
                dk, label = "undetected", "Undetected"
            elif len(dvals) == 1:
                dk = str(dvals[0])
                label = dk.split(":",1)[1] if ":" in dk else dk
            else:
                dk    = "multiple domains"
                label = " & ".join(str(d).split(":",1)[-1] for d in dvals)
        left_edge  = (left  == 1) or (left  in set(sub["qstart"].tolist()))
        right_edge = (right == qlength) or (right in set(sub["qend"].tolist()))
        ps = left  if left_edge  else left  + 1
        pe = right if right_edge else right - 1
        if ps > pe: continue
        tiles.append({"s": int(ps), "e": int(pe),
                      "dk": dk, "c": color_map.get(dk, color_map["undetected"]), "l": label})
    return tiles

protein_tiles = {}
for q, sub in resolved[resolved["query"].isin(valid_ids)].groupby("query", sort=False):
    lengths = sub["qLen"].dropna().unique()
    qlen = int(lengths[0]) if len(lengths) else int(sub["qend"].max())
    protein_tiles[q] = {"len": qlen, "tiles": build_tiles(sub, qlen)}

for pid in valid_ids:
    if pid not in protein_tiles:
        row  = valid_rbp[valid_rbp["proteinID"] == pid].iloc[0]
        plen = int(row["protein_len"])
        protein_tiles[pid] = {"len": plen,
            "tiles": [{"s":1,"e":plen,"dk":"undetected","c":color_map["undetected"],"l":"Undetected"}]}

# ── 6. Genome metadata ─────────────────────────────────────────────────────────
valid_rbp = valid_rbp.copy()
valid_rbp["Genome_name"] = valid_rbp["proteinID"].str.replace(r"(_PROTEIN_.*|_gp.*)", "", regex=True)
genome_lookup = genome_df.set_index("Genome_name").to_dict("index")

def hr_label(raw):
    if str(raw) in ("Intermediate","Broad"): return "I/B"
    if str(raw) == "Narrow": return "Narrow"
    return "Unknown"

# ── 7. Build hierarchy ─────────────────────────────────────────────────────────
PMETA = ["protein_len","RBP","PFAM_function","ECOD_function","PHROGS1_function",
         "PHROGS2_function","preds","score","ranking_score","RBPseg(YES/NO)",
         "Avg_MergeScore","Model_Quality","TF-class"]

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
        "Phage":          str(gm.get("Phage","")),
        "Genome_name":    str(gm.get("Genome_name", gname)),
        "Genome_ID":      str(gm.get("Genome_ID","")),
        "Genus":          str(gm.get("Genus","")),
        "Phage_cluster":  str(gm.get("Phage_cluster","")),
        "Genome_size_bp": str(gm.get("Genome_size_(bp)","")),
        "Morphology":     str(gm.get("Morphology","")),
        "Capsule":        str(gm.get("Capsule","")),
        "Observation":    str(gm.get("Observation","")),
        "PMID":           str(gm.get("Publication(PMID)","")),
        "Host_range":     hr_raw,
        "Halo":           str(gm.get("Halo","")),
        "HR":             hr,
    }

    ti = protein_tiles[pid]
    hierarchy[cls][cluster].append({
        "id":    pid,
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

# ── 8. Serialise hierarchy ─────────────────────────────────────────────────────
def _ckey(k):
    m2 = re.match(r"RBP(\d+)$", k)
    return int(m2.group(1)) if m2 else 9999

classes_list = []
for cls in sorted(hierarchy):
    clusters = hierarchy[cls]
    cl = [{"id":cid, "proteins":sorted(clusters[cid], key=lambda p: p["id"])}
          for cid in sorted(clusters, key=_ckey)]
    n = sum(len(c["proteins"]) for c in cl)
    classes_list.append({"name": cls, "clusters": cl, "n": n})

atlas_data = {"colorMap": color_map, "classes": classes_list}
data_json  = json.dumps(atlas_data, separators=(",",":"), ensure_ascii=False)
print(f"JSON: {len(data_json)//1024} KB  |  {len(valid_ids)} proteins  |  {len(classes_list)} classes")

# ── 9. Pseudo-domain map ───────────────────────────────────────────────────────
pd_records = []
for _, row in pd_df.iterrows():
    try:
        qTM = float(row["complex_qTMscore"])
        tTM = float(row["complex_tTMscore"])
    except (ValueError, TypeError):
        continue
    pd_records.append({
        "query":    str(row["query"]),
        "target":   str(row["target"]),
        "qTM":      round(qTM, 3),
        "tTM":      round(tTM, 3),
        "minTM":    round(min(qTM, tTM), 3),
        "tstart":   int(row["tstart"]),
        "tend":     int(row["tend"]),
        "rbpClass": str(row["RBP-class"])
    })

pd_json  = json.dumps(pd_records, separators=(",",":"), ensure_ascii=False)
pd_count = len(set(r["query"] for r in pd_records))
print(f"Pseudo-domain clusters: {pd_count}  |  Mappings: {len(pd_records)}")

# ── 9b. Sequence modularity network ────────────────────────────────────────────
sm_df = pd.read_csv(SM_MAP, sep="\t", low_memory=False)
sm_df["modular"] = sm_df["modular"].str.replace(
    "Recombination hotspots", "Putative recombination hotspots", regex=False
)

class_n = {c["name"]: c["n"] for c in classes_list}

from collections import Counter as _Counter
_edge_counts = _Counter()
for _, _row in sm_df.iterrows():
    _qc  = str(_row["query_class"])
    _tc  = str(_row["target_class"])
    _cat = str(_row["modular"])
    _edge_counts[(_qc, _tc, _cat)] += 1

sm_edges = [{"s": s, "t": t, "cat": c, "n": n}
            for (s, t, c), n in sorted(_edge_counts.items(), key=lambda x: -x[1])]

sm_classes = sorted({e["s"] for e in sm_edges} | {e["t"] for e in sm_edges})
sm_nodes   = [{"id": cls, "n": class_n.get(cls, 0)} for cls in sm_classes]

sm_data  = {"nodes": sm_nodes, "edges": sm_edges}
sm_json  = json.dumps(sm_data, separators=(",",":"), ensure_ascii=False)
print(f"SM nodes: {len(sm_nodes)}  |  SM edges (class-level): {len(sm_edges)}")

# Protein-level rows for the detail table
sm_row_records = []
for _, _row in sm_df.iterrows():
    _cat = str(_row["modular"])
    if not _cat or _cat == "nan":
        continue
    try:
        sm_row_records.append({
            "q":    str(_row["query"]),
            "t":    str(_row["target"]),
            "qc":   str(_row["query_class"]),
            "tc":   str(_row["target_class"]),
            "cat":  _cat,
            "pid":  round(float(_row["pident"]), 1),
            "al":   int(_row["alnlen"]),
            "qs":   int(_row["qstart"]), "qe": int(_row["qend"]),
            "ts":   int(_row["tstart"]), "te": int(_row["tend"]),
            "qcov": round(float(_row["qcov"]), 3),
            "tcov": round(float(_row["tcov"]), 3),
            "ev":   float(_row["evalue"]),
        })
    except (ValueError, TypeError):
        continue

sm_rows_json = json.dumps(sm_row_records, separators=(",",":"), ensure_ascii=False)
print(f"SM protein-level rows: {len(sm_row_records)}")

# ── 10. HTML template ──────────────────────────────────────────────────────────
HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Klebsiella phage RBP Atlas</title>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{display:flex;flex-direction:column;height:100vh;overflow:hidden;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;font-size:13px;background:#f4f6f9;color:#222}

/* ── Top bar ── */
#top-bar{display:flex;align-items:center;gap:12px;padding:0 20px;height:52px;min-height:52px;background:#1a2b3c;color:#fff;flex-shrink:0;z-index:50}
#top-title{font-size:16px;font-weight:700;cursor:pointer;white-space:nowrap;letter-spacing:.3px}
#top-title:hover{color:#9ac0f0}
#btn-home{display:none;align-items:center;gap:5px;background:none;border:none;color:#7a9ab8;cursor:pointer;font-size:12px;padding:4px 10px;border-radius:4px;white-space:nowrap}
#btn-home:hover{background:#253d52;color:#fff}
#btn-home.vis{display:flex}
.top-search{flex:1;max-width:360px;position:relative}
.top-search input{width:100%;padding:6px 10px;border-radius:5px;border:none;background:#253d52;color:#e0eaf4;font-size:12px;outline:none}
.top-search input::placeholder{color:#5a7a94}
.top-search input:focus{background:#2d4d66;box-shadow:0 0 0 2px #4a90d9}
#search-results{position:absolute;top:calc(100% + 3px);left:0;right:0;background:#1a2b3c;border:1px solid #2d4157;border-radius:6px;max-height:340px;overflow-y:auto;z-index:300;display:none}
.sr-section{padding:5px 14px 3px;font-size:10px;font-weight:700;color:#4a90d9;text-transform:uppercase;letter-spacing:.7px;border-top:1px solid #1e3347}
.sr-section:first-child{border-top:none}
.sr-item{padding:6px 14px;cursor:pointer;border-bottom:1px solid #1a2f42;transition:background .12s}
.sr-item:hover{background:#253d52}
.sr-pid{font-size:11px;color:#9ac0f0;display:flex;align-items:center;gap:5px}
.sr-cls{font-size:10px;color:#5a7a94;margin-top:1px}
#view-label{font-size:12px;color:#7a9ab8;margin-left:auto;white-space:nowrap}

/* ── Content wrapper ── */
#content{flex:1;display:flex;overflow:hidden}

/* ── Home ── */
#view-home{flex:1;overflow-y:auto;display:flex;flex-direction:column;align-items:center;padding:36px 20px 48px;gap:28px}
.atlas-title{font-size:34px;font-weight:700;color:#1a2b3c;text-align:center;letter-spacing:-.5px;margin-top:8px}
.atlas-sub{font-size:14px;color:#888;text-align:center;max-width:500px;line-height:1.75}
.icon-grid{display:flex;flex-wrap:wrap;gap:20px;justify-content:center}
.icon-card{background:#fff;border:2px solid #e0e6ec;border-radius:14px;padding:20px 26px 16px;text-align:center;cursor:pointer;transition:all .2s;min-width:138px;user-select:none;display:flex;flex-direction:column;align-items:center;gap:6px}
.icon-card:hover{border-color:#4a90d9;box-shadow:0 6px 20px rgba(45,90,142,.18);transform:translateY(-3px)}
.icon-card:active{transform:translateY(-1px)}
.ic-img{width:52px;height:52px;display:flex;align-items:center;justify-content:center}
.ic-num{font-size:26px;font-weight:700;color:#2d5a8e;line-height:1}
.ic-lbl{font-size:11px;color:#777;font-weight:500}

/* Donut charts section */
#home-charts{display:flex;flex-wrap:wrap;gap:32px;justify-content:center;width:100%;max-width:900px}
.chart-card{background:#fff;border:1px solid #e0e6ec;border-radius:12px;padding:20px 24px;flex:1;min-width:300px;max-width:420px}
.chart-title{font-size:13px;font-weight:700;color:#1a2b3c;margin-bottom:4px}
.chart-sub{font-size:11px;color:#aaa;margin-bottom:14px}
.donut-wrap{display:flex;align-items:center;gap:20px}
.donut-wrap canvas{flex-shrink:0}
.chart-legend{flex:1;overflow-y:auto;max-height:180px;font-size:11px}
.chart-legend::-webkit-scrollbar{width:3px}
.chart-legend::-webkit-scrollbar-thumb{background:#dde3ea;border-radius:2px}
.leg-row{display:flex;align-items:center;gap:6px;padding:2px 0;cursor:default}
.leg-dot{width:10px;height:10px;border-radius:50%;flex-shrink:0}
.leg-name{flex:1;color:#555;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.leg-pct{color:#aaa;flex-shrink:0}

/* ── Explorer ── */
#view-explorer{flex:1;overflow:hidden;display:none;flex-direction:column}
#view-hdr{display:flex;align-items:center;gap:12px;padding:10px 22px;background:#fff;border-bottom:1px solid #e0e6ec;flex-shrink:0}
#view-hdr-title{font-size:14px;font-weight:700;color:#1a2b3c;flex:1}
.hdr-btn{font-size:11px;padding:5px 12px;border-radius:5px;border:1px solid #dde3ea;background:#f8f9fa;cursor:pointer;color:#555;transition:all .15s;display:flex;align-items:center;gap:5px;white-space:nowrap}
.hdr-btn:hover{border-color:#4a90d9;color:#2d5a8e}
#explorer-content{flex:1;overflow-y:auto;padding:16px 22px}
#explorer-content::-webkit-scrollbar{width:6px}
#explorer-content::-webkit-scrollbar-thumb{background:#c8d0da;border-radius:3px}

/* Group accordion */
.grp-block{background:#fff;border:1px solid #e0e6ec;border-radius:10px;margin-bottom:10px;overflow:hidden}
.grp-head{display:flex;align-items:center;gap:10px;padding:11px 16px;cursor:pointer;user-select:none}
.grp-head:hover{background:#fafbff}
.grp-color-dot{width:13px;height:13px;border-radius:50%;flex-shrink:0}
.grp-name{font-size:13px;font-weight:600;color:#1a2b3c;flex:1;line-height:1.4}
.grp-count{font-size:11px;color:#888;white-space:nowrap}
.grp-arrow{font-size:10px;color:#bbb;transition:transform .2s;flex-shrink:0}
.grp-block.open>.grp-head>.grp-arrow{transform:rotate(90deg)}
.grp-body{display:none;border-top:1px solid #f0f3f7}
.grp-block.open>.grp-body{display:block}

/* Sub-group */
.sub-block{border-bottom:1px solid #f0f3f7}
.sub-block:last-child{border-bottom:none}
.sub-head{display:flex;align-items:center;gap:8px;padding:7px 16px 7px 24px;cursor:pointer;background:#fafbfc;user-select:none}
.sub-head:hover{background:#f0f4ff}
.sub-name{font-size:12px;font-weight:600;color:#2d5a8e;flex:1}
.sub-count{font-size:10px;color:#aaa}
.sub-arrow{font-size:9px;color:#bbb;transition:transform .2s}
.sub-block.open>.sub-head>.sub-arrow{transform:rotate(90deg)}
.sub-body{display:none}
.sub-block.open>.sub-body{display:block}

/* Protein rows */
.prot-row{display:flex;align-items:center;gap:10px;padding:5px 16px 5px 36px;cursor:pointer;border-bottom:1px solid #f8f9fb;transition:background .1s}
.prot-row:last-child{border-bottom:none}
.prot-row:hover{background:#f0f4ff}
.prot-row.selected{background:#e8f0fe;border-left:3px solid #4a90d9}
.prot-name{font-family:'SF Mono',Monaco,monospace;font-size:11px;color:#1a2b3c;flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.prot-badge{font-size:10px;font-weight:700;padding:1px 7px;border-radius:10px;flex-shrink:0}
.hr-Narrow{background:#fff0f0;color:#c0392b}
.hr-IB{background:#f0fff4;color:#27ae60}
.hr-Unknown{background:#f5f5f5;color:#999}
.prot-morph{font-size:10px;color:#aaa;flex-shrink:0;min-width:64px;text-align:right}
.prot-row-ecod{display:flex;align-items:center;gap:8px;padding:5px 16px 5px 36px;cursor:pointer;border-bottom:1px solid #f8f9fb;transition:background .1s}
.prot-row-ecod:last-child{border-bottom:none}
.prot-row-ecod:hover{background:#f0f4ff}
.prot-row-ecod.selected{background:#e8f0fe;border-left:3px solid #4a90d9}
.prot-name-ecod{font-family:'SF Mono',Monaco,monospace;font-size:10px;color:#555;min-width:180px;flex-shrink:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.domain-bar{display:flex;height:14px;border-radius:3px;overflow:hidden;border:1px solid #dde3ea;flex:1;min-width:0}
.dseg{height:100%;transition:filter .1s}
.dseg:hover{filter:brightness(.85)}
.dseg.undetected{border-right:1px solid #ececec}

/* ── Detail panel ── */
#detail{width:460px;min-width:280px;max-width:740px;background:#fff;border-left:1px solid #dde3ea;display:flex;flex-direction:column;flex-shrink:0;overflow:hidden;position:relative}
#detail.closed{width:0!important;min-width:0;border:none}
#det-resize{position:absolute;left:0;top:0;bottom:0;width:5px;cursor:ew-resize;z-index:10}
#det-resize:hover,#det-resize:active{background:rgba(74,144,217,.3)}
.det-header{padding:12px 16px 10px;border-bottom:1px solid #eef1f6;flex-shrink:0}
.det-close{float:right;background:none;border:none;cursor:pointer;font-size:18px;color:#bbb;line-height:1;padding:2px 5px;border-radius:4px}
.det-close:hover{color:#555;background:#f0f3f7}
.det-pid{font-family:'SF Mono',Monaco,monospace;font-size:12px;font-weight:700;color:#1a2b3c;padding-right:28px;margin-bottom:5px;word-break:break-all}
.det-tags{display:flex;gap:5px;flex-wrap:wrap}
.det-tag{font-size:10px;padding:2px 8px;border-radius:10px;background:#f0f3f7;color:#555}
.det-tag.cls{background:#e8f0fe;color:#2d5a8e;font-weight:600}
.det-tag.pd{background:#fff3e0;color:#e65100;font-weight:600}
.det-body{flex:1;overflow-y:auto}
.det-body::-webkit-scrollbar{width:5px}
.det-body::-webkit-scrollbar-thumb{background:#e0e6ec;border-radius:3px}
.det-section{padding:12px 16px;border-bottom:1px solid #eef1f6}
.det-section h3{font-size:11px;font-weight:700;color:#999;text-transform:uppercase;letter-spacing:.6px;margin-bottom:9px}

/* NGL */
#ngl-viewport{height:250px;background:#f8f9fa;border-radius:6px;overflow:hidden;border:1px solid #e0e6ec;position:relative;resize:vertical;min-height:120px}
.ngl-loading{position:absolute;inset:0;display:flex;align-items:center;justify-content:center;font-size:12px;color:#bbb}
.ngl-controls{display:flex;gap:6px;margin-top:8px;flex-wrap:wrap;align-items:center}
.ngl-btn{font-size:11px;padding:3px 10px;border-radius:4px;border:1px solid #dde3ea;background:#f8f9fa;cursor:pointer;color:#555;transition:all .15s;white-space:nowrap}
.ngl-btn:hover{border-color:#4a90d9;color:#2d5a8e}
.ngl-btn.active{background:#2d5a8e;color:#fff;border-color:#2d5a8e}
.ngl-no-struct{display:flex;align-items:center;justify-content:center;height:100%;color:#ccc;font-size:12px}
a.ngl-btn{text-decoration:none;display:inline-flex;align-items:center;gap:4px}

/* Fullscreen overlay */
#ngl-fullscreen{display:none;position:fixed;inset:0;z-index:9000;background:#fff;flex-direction:column}
#ngl-fullscreen.active{display:flex}
#ngl-fs-bar{display:flex;align-items:center;gap:10px;padding:8px 14px;background:#1a2b3c;color:#fff;flex-shrink:0}
#ngl-fs-bar span{flex:1;font-size:13px;font-weight:600}
#ngl-fs-viewport{flex:1}
#btn-fs-close{background:none;border:none;color:#aaa;cursor:pointer;font-size:20px;padding:2px 6px;border-radius:4px}
#btn-fs-close:hover{color:#fff;background:#253d52}

/* Sequence strip */
#sec-seq{padding:10px 16px 12px}
#seq-strip-wrap{position:relative}
#seq-canvas{display:block;width:100%;border-radius:4px;cursor:crosshair}
#seq-pos-label{font-size:10px;color:#888;margin-top:3px;height:14px;font-family:monospace}

/* Domain architecture */
.domain-bar-large{display:flex;height:24px;border-radius:4px;overflow:hidden;border:1px solid #dde3ea}
.bar-ruler{display:flex;justify-content:space-between;font-size:10px;color:#bbb;margin-top:2px;font-family:monospace}
.domain-legend{display:flex;flex-wrap:wrap;gap:7px;margin-top:8px}
.leg-item{display:flex;align-items:center;gap:4px;font-size:11px;color:#555}
.leg-swatch{width:13px;height:13px;border-radius:2px;border:1px solid rgba(0,0,0,.1);flex-shrink:0}

/* PD region bar */
.pd-bar-wrap{position:relative;height:26px;border-radius:4px;background:#f0f3f7;border:1px solid #dde3ea;overflow:hidden;margin-bottom:3px}
.pd-seg{position:absolute;top:0;bottom:0;opacity:.82}

/* Metadata tables */
.tabs{display:flex;border-bottom:2px solid #eef1f6;margin-bottom:10px}
.tab-btn{padding:6px 14px;font-size:12px;font-weight:600;color:#aaa;border:none;background:none;cursor:pointer;border-bottom:2px solid transparent;margin-bottom:-2px;transition:all .15s}
.tab-btn:hover{color:#2d5a8e}
.tab-btn.active{color:#2d5a8e;border-bottom-color:#2d5a8e}
.tab-content{display:none}
.tab-content.active{display:block}
.meta-table{width:100%;border-collapse:collapse;font-size:11.5px}
.meta-table tr:not(:last-child) td{border-bottom:1px solid #f4f6f9}
.meta-table td{padding:4px 6px;vertical-align:top}
.meta-table td:first-child{color:#999;font-weight:600;white-space:nowrap;width:130px}
.meta-table td:last-child{color:#333;word-break:break-word}
.pmid-link{color:#4a90d9;text-decoration:none}
.pmid-link:hover{text-decoration:underline}
.tm-badge{display:inline-block;padding:1px 6px;border-radius:8px;font-size:10px;font-weight:700}
.tm-h{background:#e8f5e9;color:#2e7d32}
.tm-m{background:#fff8e1;color:#f57f17}
.tm-l{background:#fce4ec;color:#c62828}

/* Tooltip */
#tooltip{position:fixed;pointer-events:none;background:rgba(15,25,40,.88);color:#e8f0f8;font-size:11px;padding:5px 9px;border-radius:5px;max-width:260px;line-height:1.5;z-index:9999;display:none}

/* ── Modularity network view ── */
#view-modularity{flex:1;overflow:hidden;display:none;flex-direction:column}
#net-hdr{display:flex;align-items:center;gap:10px;padding:10px 22px;background:#fff;border-bottom:1px solid #e0e6ec;flex-shrink:0;flex-wrap:wrap}
.net-legend{display:flex;flex-wrap:wrap;gap:7px;align-items:center;flex:1}
.cat-chip{display:inline-flex;align-items:center;gap:5px;font-size:11px;padding:3px 9px;border-radius:10px;background:#f4f6f9;color:#444;white-space:nowrap;border:1px solid #e0e6ec}
.cat-swatch{width:10px;height:10px;border-radius:50%;flex-shrink:0}
#cy-container{flex:1;min-height:0;background:#fafbfc}
#net-table-panel{flex-shrink:0;max-height:260px;overflow-y:auto;border-top:2px solid #e0e6ec;background:#fff;padding:12px 22px;display:none}
#net-table-panel::-webkit-scrollbar{width:5px}
#net-table-panel::-webkit-scrollbar-thumb{background:#e0e6ec;border-radius:3px}
.net-table{width:100%;border-collapse:collapse;font-size:11.5px}
.net-table th{font-size:10px;font-weight:700;color:#999;text-transform:uppercase;letter-spacing:.6px;border-bottom:2px solid #eef1f6;padding:4px 8px;text-align:left}
.net-table td{padding:5px 8px;border-bottom:1px solid #f4f6f9;vertical-align:middle}
</style>
</head>
<body>

<header id="top-bar">
  <span id="top-title" onclick="goHome()">Klebsiella phage RBP Atlas</span>
  <button id="btn-home" onclick="goHome()">&#8592; Home</button>
  <div class="top-search">
    <input id="search-input" type="text" placeholder="Search protein, phage, capsule, genus, cluster…" autocomplete="off">
    <div id="search-results"></div>
  </div>
  <span id="view-label"></span>
</header>

<div id="content">

  <!-- ── Home ── -->
  <div id="view-home">
    <div class="atlas-title">Klebsiella phage RBP Atlas</div>
    <div class="atlas-sub">A structural atlas of receptor-binding proteins from <em>Klebsiella</em>-infecting bacteriophages. Explore structures, domain architectures, pseudo-domain diversity, and sequence modularity across RBP classes.</div>
    <div class="icon-grid" id="icon-grid"></div>
    <div id="home-charts"></div>
  </div>

  <!-- ── Explorer ── -->
  <div id="view-explorer">
    <div id="view-hdr">
      <span id="view-hdr-title"></span>
      <button class="hdr-btn" onclick="exportTSV()" title="Download visible proteins as TSV">
        &#8659; Export TSV
      </button>
    </div>
    <div id="explorer-content"></div>
  </div>

  <!-- ── Detail panel ── -->
  <aside id="detail" class="closed">
    <div id="det-resize"></div>
    <div class="det-header">
      <button class="det-close" id="close-detail">&#10005;</button>
      <div class="det-pid" id="det-pid"></div>
      <div class="det-tags" id="det-tags"></div>
    </div>
    <div class="det-body">

      <!-- 3D structure -->
      <div class="det-section">
        <h3 id="det-struct-title">3D Structure</h3>
        <div id="ngl-viewport"><div class="ngl-loading">Loading…</div></div>
        <div class="ngl-controls">
          <button class="ngl-btn" id="btn-spin">&#8635; Spin</button>
          <button class="ngl-btn" id="btn-reset">&#8853; Reset</button>
          <button class="ngl-btn" id="btn-surf">&#9724; Surface</button>
          <button class="ngl-btn" id="btn-fs" title="Fullscreen">&#10697;</button>
          <a  class="ngl-btn" id="btn-dl" download title="Download structure file">&#8659; .pdb.gz</a>
        </div>
      </div>

      <!-- Sequence strip -->
      <div class="det-section" id="sec-seq" style="display:none">
        <h3>Sequence Coverage</h3>
        <div id="seq-strip-wrap">
          <canvas id="seq-canvas" height="36"></canvas>
          <div id="seq-pos-label"></div>
        </div>
      </div>

      <!-- Domain architecture — ECOD view -->
      <div class="det-section" id="sec-arch" style="display:none">
        <h3>Domain Architecture</h3>
        <div id="det-bar"></div>
        <div id="det-ruler" class="bar-ruler"></div>
        <div id="det-legend" class="domain-legend"></div>
      </div>

      <!-- PD region — PD view -->
      <div class="det-section" id="sec-pd" style="display:none">
        <h3>Pseudo-domain Region</h3>
        <div id="det-pd-bar-wrap"></div>
        <div id="det-pd-info"></div>
      </div>

      <!-- Metadata tabs -->
      <div class="det-section">
        <div class="tabs">
          <button class="tab-btn active" data-tab="protein">Protein</button>
          <button class="tab-btn" data-tab="genome">Genome</button>
        </div>
        <div id="tab-protein" class="tab-content active"></div>
        <div id="tab-genome"  class="tab-content"></div>
      </div>

    </div><!-- det-body -->
  </aside>

  <!-- ── Modularity network ── -->
  <div id="view-modularity">
    <div id="net-hdr">
      <div class="net-legend" id="net-legend"></div>
      <button class="hdr-btn" id="btn-net-export" onclick="exportNetRowsTSV()">&#8659; Export TSV</button>
    </div>
    <div id="cy-container"></div>
    <div id="net-table-panel"></div>
  </div>

</div><!-- content -->

<!-- Fullscreen NGL overlay -->
<div id="ngl-fullscreen">
  <div id="ngl-fs-bar">
    <span id="ngl-fs-title"></span>
    <button id="btn-fs-close" title="Exit fullscreen">&#10005;</button>
  </div>
  <div id="ngl-fs-viewport"></div>
</div>

<div id="tooltip"></div>

<script src="https://cdn.jsdelivr.net/npm/ngl@0.10.4/dist/ngl.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.3/dist/chart.umd.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/cytoscape@3.28.1/dist/cytoscape.min.js"></script>
<script>
const D       = __ATLAS_DATA__;
const PD      = __PD_DATA__;
const SM      = __SM_DATA__;
const SM_ROWS = __SM_ROWS__;

const CAT_COLORS = {
  "Cterminal variation":             "#1f77b4",
  "Nterminal sharing":               "#ff7f0e",
  "Putative recombination hotspots": "#2ca02c",
  "Cterminal sharing":               "#d62728",
  "Conserved regions":               "#9467bd"
};

// ── Lookups ────────────────────────────────────────────────────────────────────
const classMap   = new Map(D.classes.map(c => [c.name, c]));
const proteinMap = new Map();
D.classes.forEach(c => c.clusters.forEach(cl =>
  cl.proteins.forEach(p => proteinMap.set(p.id, p))));

const pdByTarget = new Map();
const pdByQuery  = new Map();
PD.forEach(r => {
  if (!pdByTarget.has(r.target)) pdByTarget.set(r.target, []);
  pdByTarget.get(r.target).push(r);
  if (!pdByQuery.has(r.query))   pdByQuery.set(r.query, []);
  pdByQuery.get(r.query).push(r);
});

const PD_IDS = [...new Set(PD.map(r => r.query))].sort();
const PD_COLOR = {};
PD_IDS.forEach((id, i) => {
  PD_COLOR[id] = `hsl(${Math.round(i * 360 / PD_IDS.length)},70%,50%)`;
});

// ── State ──────────────────────────────────────────────────────────────────────
const S = {
  view: "home",
  nglStage: null, nglComp: null,
  spinning: false, surfaceRep: null, hoverRep: null,
  detailMode: null, activeProt: null, activePd: null,
  currentP: null,   // protein object currently in detail
  fsStage: null,    // fullscreen NGL stage
  cy: null,         // Cytoscape instance
  netFilter: null,  // current network table filter
};

// ── Init ──────────────────────────────────────────────────────────────────────
document.addEventListener("DOMContentLoaded", () => {
  const total         = D.classes.reduce((a,c) => a + c.n, 0);
  const totalClusters = D.classes.reduce((a,c) => a + c.clusters.length, 0);
  const ecodSet       = new Set();
  D.classes.forEach(c => c.clusters.forEach(cl => cl.proteins.forEach(p =>
    p.tiles.forEach(t => { if (t.dk !== "undetected" && t.dk !== "multiple domains") ecodSet.add(t.dk); })
  )));

  const rbpSvg = `<svg width="52" height="52" viewBox="0 0 52 52" fill="none" xmlns="http://www.w3.org/2000/svg"><path d="M10 42 C10 34 18 30 18 26 C18 22 10 18 12 10" stroke="#4a90d9" stroke-width="4.5" stroke-linecap="round"/><path d="M20 42 C20 34 28 30 28 26 C28 22 20 18 22 10" stroke="#2d5a8e" stroke-width="4.5" stroke-linecap="round"/><path d="M32 16 L46 16 L43 12 L48 16 L43 20 Z" fill="#e65100"/><path d="M32 26 L46 26 L43 22 L48 26 L43 30 Z" fill="#4a90d9"/><path d="M32 36 L46 36 L43 32 L48 36 L43 40 Z" fill="#2d5a8e"/><path d="M28 10 Q30 14 32 16" stroke="#7a9ab8" stroke-width="1.8" fill="none"/><path d="M28 26 Q30 26 32 26" stroke="#7a9ab8" stroke-width="1.8" fill="none"/><path d="M28 42 Q30 39 32 36" stroke="#7a9ab8" stroke-width="1.8" fill="none"/></svg>`;
  const classSvg = `<svg width="52" height="52" viewBox="0 0 52 52" fill="none"><rect x="6" y="8" width="40" height="9" rx="4" fill="#2d5a8e"/><rect x="6" y="22" width="28" height="9" rx="4" fill="#4a90d9"/><rect x="6" y="36" width="34" height="9" rx="4" fill="#7ab4e8"/></svg>`;
  const clusterSvg = `<svg width="52" height="52" viewBox="0 0 52 52" fill="none"><circle cx="26" cy="20" r="7" fill="#2d5a8e"/><circle cx="14" cy="36" r="5" fill="#4a90d9"/><circle cx="26" cy="36" r="5" fill="#4a90d9"/><circle cx="38" cy="36" r="5" fill="#4a90d9"/><line x1="26" y1="27" x2="14" y2="31" stroke="#b0c8e8" stroke-width="1.8"/><line x1="26" y1="27" x2="26" y2="31" stroke="#b0c8e8" stroke-width="1.8"/><line x1="26" y1="27" x2="38" y2="31" stroke="#b0c8e8" stroke-width="1.8"/></svg>`;
  const ecodSvg = `<svg width="52" height="52" viewBox="0 0 52 52" fill="none"><rect x="6" y="22" width="40" height="10" rx="3" fill="#e0e6ec"/><rect x="6" y="22" width="12" height="10" rx="3" fill="#cc0099"/><rect x="18" y="22" width="10" height="10" fill="#0047ab"/><rect x="28" y="22" width="8" height="10" fill="#e68c00"/><rect x="36" y="22" width="10" height="10" rx="3" fill="#00998c"/><path d="M26 10 L26 22 M26 32 L26 42" stroke="#bbb" stroke-width="1.5" stroke-dasharray="3,2"/></svg>`;
  const pdSvg  = `<svg width="52" height="52" viewBox="0 0 52 52" fill="none"><ellipse cx="26" cy="26" rx="20" ry="20" fill="none" stroke="#e0e6ec" stroke-width="2"/><path d="M26 6 A20 20 0 0 1 46 26" stroke="#e65100" stroke-width="5" stroke-linecap="round"/><path d="M46 26 A20 20 0 0 1 26 46" stroke="#4a90d9" stroke-width="5" stroke-linecap="round"/><path d="M26 46 A20 20 0 0 1 6 26" stroke="#2d5a8e" stroke-width="5" stroke-linecap="round"/><path d="M6 26 A20 20 0 0 1 26 6" stroke="#27ae60" stroke-width="5" stroke-linecap="round"/><circle cx="26" cy="26" r="5" fill="#1a2b3c"/></svg>`;
  const modSvg = `<svg width="52" height="52" viewBox="0 0 52 52" fill="none"><line x1="26" y1="26" x2="10" y2="13" stroke="#1f77b4" stroke-width="2.5"/><line x1="26" y1="26" x2="42" y2="13" stroke="#ff7f0e" stroke-width="2.5"/><line x1="26" y1="26" x2="10" y2="42" stroke="#2ca02c" stroke-width="2.5"/><line x1="26" y1="26" x2="42" y2="42" stroke="#d62728" stroke-width="2.5"/><line x1="10" y1="13" x2="42" y2="13" stroke="#9467bd" stroke-width="1.8"/><line x1="10" y1="13" x2="10" y2="42" stroke="#9467bd" stroke-width="1.8"/><circle cx="26" cy="26" r="7" fill="#1a2b3c"/><circle cx="10" cy="13" r="5" fill="#4a90d9"/><circle cx="42" cy="13" r="5" fill="#4a90d9"/><circle cx="10" cy="42" r="5" fill="#4a90d9"/><circle cx="42" cy="42" r="5" fill="#4a90d9"/></svg>`;

  [
    { svg: rbpSvg,     num: total,            lbl: "RBPs",                   view: "rbp" },
    { svg: classSvg,   num: D.classes.length, lbl: "RBP Classes",            view: "classes" },
    { svg: clusterSvg, num: totalClusters,    lbl: "RBP Clusters",           view: "clusters" },
    { svg: ecodSvg,    num: ecodSet.size,     lbl: "ECOD Domains",           view: "ecod" },
    { svg: pdSvg,      num: PD_IDS.length,    lbl: "Pseudo-domain Clusters", view: "pd" },
    { svg: modSvg,     num: SM.edges.length,  lbl: "Sequence Modularity",    view: "modularity" },
  ].forEach(d => {
    const el = document.createElement("div");
    el.className = "icon-card";
    el.innerHTML = `<div class="ic-img">${d.svg}</div><div class="ic-num">${d.num}</div><div class="ic-lbl">${d.lbl}</div>`;
    el.onclick = () => showView(d.view);
    document.getElementById("icon-grid").appendChild(el);
  });

  buildHomeCharts();
  setupDetailButtons();
  setupSearch();
  setupTooltip();
  setupResize();
  setupFullscreen();
});

// ── Donut charts ──────────────────────────────────────────────────────────────
function buildHomeCharts() {
  const ecodRes = new Map();
  D.classes.forEach(c => c.clusters.forEach(cl => cl.proteins.forEach(p => {
    if (!p.pdb) return;
    p.tiles.forEach(t => {
      const n = t.e - t.s + 1;
      ecodRes.set(t.dk, (ecodRes.get(t.dk)||0) + n);
    });
  })));
  const ecodEntries = [...ecodRes.entries()].sort((a,b) => b[1]-a[1]);
  const ecodTotal   = ecodEntries.reduce((s,[,v])=>s+v,0);

  const wrap = document.getElementById("home-charts");
  wrap.innerHTML = `
    <div class="chart-card">
      <div class="chart-title">ECOD Domain Coverage</div>
      <div class="chart-sub">Fraction of structural residues annotated by ECOD</div>
      <div class="donut-wrap">
        <canvas id="ecod-chart" width="160" height="160"></canvas>
        <div class="chart-legend" id="ecod-legend"></div>
      </div>
    </div>`;

  const ecodColors = ecodEntries.map(([dk]) =>
    dk === "undetected" ? "#e8e8e8" : (D.colorMap[dk] || "#aaa"));
  const ecodLabels = ecodEntries.map(([dk]) =>
    dk === "undetected" ? "Undetected" : dk.split(":")[1] || dk);

  new Chart(document.getElementById("ecod-chart"), {
    type: "doughnut",
    data: { labels: ecodLabels, datasets: [{ data: ecodEntries.map(([,v])=>v), backgroundColor: ecodColors, borderWidth: 1, borderColor: "#fff" }] },
    options: { plugins: { legend: { display: false }, tooltip: { callbacks: {
      label: (ctx) => ` ${ctx.label}: ${((ctx.raw/ecodTotal)*100).toFixed(1)}% (${ctx.raw.toLocaleString()} aa)`
    }}}, cutout: "60%", animation: { duration: 600 } }
  });
  const ecodLeg = document.getElementById("ecod-legend");
  ecodEntries.forEach(([dk, v], i) => {
    const pct = ((v/ecodTotal)*100).toFixed(1);
    const lbl = dk === "undetected" ? "Undetected" : dk.split(":")[1] || dk;
    ecodLeg.innerHTML += `<div class="leg-row"><div class="leg-dot" style="background:${ecodColors[i]}"></div>
      <span class="leg-name" title="${esc(lbl)}">${esc(lbl)}</span>
      <span class="leg-pct">${pct}%</span></div>`;
  });
}

// ── View switching ─────────────────────────────────────────────────────────────
const VIEW_LABELS = {
  rbp:"RBP Structures", classes:"RBP Classes", clusters:"RBP Clusters",
  ecod:"ECOD Domains", pd:"Pseudo-domain Clusters", modularity:"Sequence Modularity"
};

function setView(name) {
  document.getElementById("view-home").style.display       = name === "home"       ? "" : "none";
  document.getElementById("view-explorer").style.display   = (name !== "home" && name !== "modularity") ? "flex" : "none";
  document.getElementById("view-modularity").style.display = name === "modularity" ? "flex" : "none";
  document.getElementById("btn-home").classList.toggle("vis", name !== "home");
  document.getElementById("view-label").textContent = VIEW_LABELS[name] || "";
  if (name !== "home" && name !== "modularity")
    document.getElementById("view-hdr-title").textContent = VIEW_LABELS[name] || "";
  S.view = name;
}

function goHome() {
  setView("home");
  closeDetail();
  document.getElementById("explorer-content").innerHTML = "";
  if (S.cy) { S.cy.destroy(); S.cy = null; }
}

function showView(name) {
  setView(name);
  closeDetail();
  if (name === "modularity") { buildModularityView(); return; }
  const ec = document.getElementById("explorer-content");
  if      (name === "rbp")      ec.innerHTML = buildRbpView();
  else if (name === "classes")  ec.innerHTML = buildClassesView();
  else if (name === "clusters") ec.innerHTML = buildClustersView();
  else if (name === "ecod")     ec.innerHTML = buildEcodView();
  else if (name === "pd")       ec.innerHTML = buildPdView();
}

// ── TSV Export ────────────────────────────────────────────────────────────────
function exportTSV() {
  const rows = [];
  const hdr  = ["proteinID","RBP_class","RBP_cluster","morphology","host_range",
                 "genus","phage","protein_length","model_quality","PDB_available"];
  if (S.view === "pd") hdr.push("PD_cluster","tstart","tend","qTMscore","tTMscore","minTMscore");
  rows.push(hdr.join("\t"));

  const seen = new Set();
  const addProt = (p, extra=[]) => {
    const key = p.id + (extra.join("|")||"");
    if (seen.has(key)) return; seen.add(key);
    rows.push([
      p.id, classOfProtein(p.id)||"", p.clid||"", p.morph||"", p.hr||"",
      p.gm.Genus||"", p.gm.Phage||"", p.plen||"",
      p.pm.Model_Quality||"", p.pdb?"yes":"no",
      ...extra
    ].join("\t"));
  };

  if (S.view === "rbp" || S.view === "classes" || S.view === "clusters") {
    D.classes.forEach(c => c.clusters.forEach(cl => cl.proteins.forEach(p => addProt(p))));
  } else if (S.view === "ecod") {
    const seen2 = new Set();
    D.classes.forEach(c => c.clusters.forEach(cl => cl.proteins.forEach(p => {
      if (seen2.has(p.id)) return; seen2.add(p.id);
      if (p.tiles.some(t => t.dk !== "undetected" && t.dk !== "multiple domains")) addProt(p);
    })));
  } else if (S.view === "pd") {
    PD.forEach(r => {
      const p = proteinMap.get(r.target);
      if (p) addProt(p, [r.query, r.tstart, r.tend, r.qTM, r.tTM, r.minTM]);
    });
  }

  const blob = new Blob([rows.join("\n")], { type: "text/tab-separated-values" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `rbp-atlas-${S.view}-${new Date().toISOString().slice(0,10)}.tsv`;
  a.click(); URL.revokeObjectURL(a.href);
}

// ── View builders ──────────────────────────────────────────────────────────────
function buildRbpView() {
  const byMorph = new Map();
  D.classes.forEach(c => c.clusters.forEach(cl => cl.proteins.forEach(p => {
    const m = p.morph || "Unknown";
    if (!byMorph.has(m)) byMorph.set(m, []);
    byMorph.get(m).push(p);
  })));
  const MC = { Myovirus:"#2d5a8e", Siphovirus:"#27ae60", Podovirus:"#e65100", Unknown:"#888" };
  return [...byMorph.entries()].sort((a,b)=>b[1].length-a[1].length)
    .map(([m,ps]) => grpBlock(`morph-${esc(m)}`,
      `<span class="grp-color-dot" style="background:${MC[m]||"#888"}"></span>
       <span class="grp-name">${esc(m)}</span>
       <span class="grp-count">${ps.length} RBP${ps.length>1?"s":""}</span>`,
      ps.map(p=>protRowBasic(p,"basic")).join(""), true)).join("");
}

function buildClassesView() {
  return D.classes.map(c => {
    const all = c.clusters.flatMap(cl=>cl.proteins);
    return grpBlock(`cls-${esc(c.name)}`,
      `<span class="grp-name">${esc(c.name)}</span>
       <span class="grp-count">${all.length} protein${all.length>1?"s":""}</span>`,
      all.map(p=>protRowBasic(p,"basic")).join(""));
  }).join("");
}

function buildClustersView() {
  return D.classes.flatMap(c => c.clusters.map(cl =>
    grpBlock(`cl-${esc(cl.id)}`,
      `<span class="grp-name">${esc(cl.id)}</span>
       <span class="grp-count">${cl.proteins.length} · <em style="color:#aaa">${esc(c.name)}</em></span>`,
      cl.proteins.map(p=>protRowBasic(p,"basic")).join(""))
  )).join("");
}

function buildEcodView() {
  const ecodMap = new Map();
  D.classes.forEach(c => c.clusters.forEach(cl => cl.proteins.forEach(p => {
    const seen = new Set();
    p.tiles.forEach(t => {
      if (t.dk==="undetected"||t.dk==="multiple domains") return;
      if (seen.has(t.dk)) return; seen.add(t.dk);
      if (!ecodMap.has(t.dk)) ecodMap.set(t.dk,{label:t.l,color:t.c,byClass:new Map()});
      const g = ecodMap.get(t.dk);
      if (!g.byClass.has(c.name)) g.byClass.set(c.name,[]);
      g.byClass.get(c.name).push(p);
    });
  })));
  return [...ecodMap.entries()]
    .sort((a,b)=>{
      const ta=[...a[1].byClass.values()].reduce((s,v)=>s+v.length,0);
      const tb=[...b[1].byClass.values()].reduce((s,v)=>s+v.length,0);
      return tb-ta;
    })
    .map(([dk,info]) => {
      const total = [...info.byClass.values()].reduce((s,v)=>s+v.length,0);
      let sub=""; info.byClass.forEach((ps,cls) => {
        sub += subBlock(`ecod-${esc(dk)}-${esc(cls)}`,
          `<span class="sub-name">${esc(cls)}</span><span class="sub-count">${ps.length}</span>`,
          ps.map(p=>protRowEcod(p)).join(""));
      });
      return grpBlock(`ecod-${esc(dk)}`,
        `<span class="grp-color-dot" style="background:${info.color}"></span>
         <span class="grp-name">${esc(info.label)}</span>
         <span class="grp-count">${total}</span>`, sub);
    }).join("");
}

function buildPdView() {
  return PD_IDS.map(pdId => {
    const mappings = pdByQuery.get(pdId)||[];
    const byClass  = new Map();
    mappings.forEach(r => {
      const cls = classOfProtein(r.target)||r.rbpClass;
      if (!byClass.has(cls)) byClass.set(cls,[]);
      byClass.get(cls).push(r);
    });
    const col = PD_COLOR[pdId];
    let sub=""; byClass.forEach((ms,cls) => {
      sub += subBlock(`pd-${esc(pdId)}-${esc(cls)}`,
        `<span class="sub-name">${esc(cls)}</span><span class="sub-count">${ms.length}</span>`,
        ms.map(r=>protRowPd(r,pdId)).join(""));
    });
    return grpBlock(`pdgrp-${esc(pdId)}`,
      `<span class="grp-color-dot" style="background:${col}"></span>
       <span class="grp-name">${esc(pdId)}</span>
       <span class="grp-count">${mappings.length} structure${mappings.length>1?"s":""}</span>`, sub);
  }).join("");
}

// ── Row helpers ────────────────────────────────────────────────────────────────
function protRowBasic(p,mode) {
  const hrC = "hr-"+p.hr.replace("/","");
  const sel = S.activeProt===p.id&&S.detailMode===mode?" selected":"";
  return `<div class="prot-row${sel}" id="pr-${esc(p.id)}" onclick="openDetail('${esc(p.id)}','${mode}')">
    <span class="prot-name">${esc(p.id)}</span>
    <span class="prot-badge ${hrC}">${esc(p.hr)}</span>
    <span class="prot-morph">${esc(p.morph||"")}</span>
  </div>`;
}
function protRowEcod(p) {
  const sel = S.activeProt===p.id&&S.detailMode==="ecod"?" selected":"";
  return `<div class="prot-row-ecod${sel}" id="pr-ecod-${esc(p.id)}" onclick="openDetail('${esc(p.id)}','ecod')">
    <span class="prot-name-ecod" title="${esc(p.id)}">${esc(p.id)}</span>
    ${domainBarHTML(p,14)}
  </div>`;
}
function protRowPd(r,pdId) {
  const sel = S.activeProt===r.target&&S.activePd===pdId?" selected":"";
  return `<div class="prot-row${sel}" id="pr-pd-${esc(pdId)}-${esc(r.target)}" onclick="openDetail('${esc(r.target)}','pd','${esc(pdId)}')">
    <span class="prot-name">${esc(r.target)}</span>
    <span class="tm-badge ${tmCls(r.minTM)}" style="margin-left:auto">${r.minTM}</span>
  </div>`;
}
function tmCls(v){ return v>=0.5?"tm-h":v>=0.3?"tm-m":"tm-l"; }

// ── Accordion ─────────────────────────────────────────────────────────────────
function domainBarHTML(p,h) {
  return `<div class="domain-bar" style="height:${h}px">${p.tiles.map(t=>{
    const b=t.dk==="undetected"?"border-right:1px solid #ececec;":"";
    return `<div class="dseg${t.dk==="undetected"?" undetected":""}" style="flex:${t.e-t.s+1};background:${t.c};${b}height:${h}px" data-tip="${esc(t.l)} (${t.s}–${t.e})"></div>`;
  }).join("")}</div>`;
}
function grpBlock(id,hd,body,open=false){
  return `<div class="grp-block${open?" open":""}" id="gb-${esc(id)}">
    <div class="grp-head" onclick="toggleGrp('${esc(id)}')">${hd}<span class="grp-arrow">&#9658;</span></div>
    <div class="grp-body">${body}</div></div>`;
}
function subBlock(id,hd,body){
  return `<div class="sub-block" id="sb-${esc(id)}">
    <div class="sub-head" onclick="toggleSub('${esc(id)}')">${hd}<span class="sub-arrow">&#9658;</span></div>
    <div class="sub-body">${body}</div></div>`;
}
function toggleGrp(id){ document.getElementById("gb-"+id)?.classList.toggle("open"); }
function toggleSub(id){ document.getElementById("sb-"+id)?.classList.toggle("open"); }

// ── Detail panel ──────────────────────────────────────────────────────────────
function openDetail(protId, mode, pdId=null) {
  const p = proteinMap.get(protId);
  if (!p) return;
  S.activeProt=protId; S.detailMode=mode; S.activePd=pdId; S.currentP=p;

  document.querySelectorAll(".prot-row,.prot-row-ecod").forEach(el=>el.classList.remove("selected"));
  const rowId = mode==="ecod"?`pr-ecod-${protId}`:mode==="pd"?`pr-pd-${pdId}-${protId}`:`pr-${protId}`;
  document.getElementById(rowId)?.classList.add("selected");

  document.getElementById("detail").classList.remove("closed");
  document.getElementById("det-pid").textContent = p.id;

  const tags = document.getElementById("det-tags");
  tags.innerHTML="";
  const addTag=(t,c="")=>t&&(tags.innerHTML+=`<span class="det-tag ${c}">${esc(String(t))}</span>`);
  addTag(classOfProtein(protId),"cls");
  if (pdId) addTag(pdId,"pd");
  addTag(p.hr); addTag(p.morph); addTag(p.plen+" aa");

  document.getElementById("sec-arch").style.display = mode==="ecod"?"":"none";
  document.getElementById("sec-pd").style.display   = mode==="pd" ?"":"none";
  document.getElementById("sec-seq").style.display  = p.pdb       ?"":"none";

  // ECOD architecture
  if (mode==="ecod") {
    document.getElementById("det-bar").innerHTML =
      `<div class="domain-bar-large">${p.tiles.map(t=>
        `<div class="dseg${t.dk==="undetected"?" undetected":""}" style="flex:${t.e-t.s+1};background:${t.c};height:24px" data-tip="${esc(t.l)} (${t.s}–${t.e})"></div>`).join("")}</div>`;
    document.getElementById("det-ruler").innerHTML=`<span>1</span><span>${p.tlen}</span>`;
    const seen=new Map(); p.tiles.forEach(t=>{if(!seen.has(t.dk))seen.set(t.dk,t);});
    document.getElementById("det-legend").innerHTML=[...seen.values()].filter(t=>t.dk!=="undetected")
      .map(t=>`<div class="leg-item"><div class="leg-swatch" style="background:${t.c}"></div><span>${esc(t.l)}</span></div>`).join("");
  }

  // PD info
  if (mode==="pd"&&pdId) {
    const map=(pdByQuery.get(pdId)||[]).find(r=>r.target===protId);
    if (map) {
      const len=p.tlen||p.plen||1;
      const col=PD_COLOR[pdId]||"#f60";
      const l=((map.tstart-1)/len*100).toFixed(2), w=((map.tend-map.tstart+1)/len*100).toFixed(2);
      document.getElementById("det-pd-bar-wrap").innerHTML=
        `<div class="pd-bar-wrap" data-tip="${esc(pdId)}: ${map.tstart}–${map.tend}">
           <div class="pd-seg" style="left:${l}%;width:${w}%;background:${col}"></div>
         </div><div class="bar-ruler"><span>1</span><span>${len}</span></div>`;
      document.getElementById("det-pd-info").innerHTML=`
        <table class="meta-table" style="margin-top:8px">
          <tr><td>PD cluster</td><td><span style="display:inline-flex;align-items:center;gap:5px">
            <span style="width:10px;height:10px;border-radius:50%;background:${col};display:inline-block"></span>${esc(pdId)}</span></td></tr>
          <tr><td>Region</td><td>${map.tstart} – ${map.tend}</td></tr>
          <tr><td>qTMscore</td><td>${map.qTM}</td></tr>
          <tr><td>tTMscore</td><td>${map.tTM}</td></tr>
          <tr><td>min TMscore</td><td><span class="tm-badge ${tmCls(map.minTM)}">${map.minTM}</span></td></tr>
          <tr><td>RBP class</td><td>${esc(map.rbpClass)}</td></tr>
        </table>`;
    }
  }

  // Download link
  const dlBtn = document.getElementById("btn-dl");
  if (p.pdb) {
    dlBtn.href = p.pdb;
    dlBtn.setAttribute("download", p.id+".pdb.gz");
    dlBtn.style.display="";
  } else {
    dlBtn.style.display="none";
  }

  // Metadata
  const PM={protein_len:"Length (aa)",RBP:"RBP prediction",PFAM_function:"Pfam",ECOD_function:"ECOD",
    PHROGS1_function:"PHROGS1",PHROGS2_function:"PHROGS2",preds:"Confidence",score:"Score",
    ranking_score:"Ranking score",RBPseg:"RBPseg",Avg_MergeScore:"Avg merge score",
    Model_Quality:"Model quality","TF-class":"TF class"};
  const GM={Phage:"Phage name",Genome_name:"Genome",Genome_ID:"Genome ID",Genus:"Genus",
    Phage_cluster:"Phage cluster",Genome_size_bp:"Genome size (bp)",Morphology:"Morphology",
    Capsule:"Capsule",Observation:"Observation",PMID:"PubMed ID",Host_range:"Host range",
    Halo:"Halo",HR:"HR group"};
  document.getElementById("tab-protein").innerHTML=makeTable(p.pm,PM);
  document.getElementById("tab-genome").innerHTML =makeTable(p.gm,GM);

  const note=mode==="ecod"?" (ECOD colors)":mode==="pd"?` (${pdId} highlighted)`:" (B-factor)";
  document.getElementById("det-struct-title").textContent="3D Structure"+note;

  initNGL(p,mode,pdId);
}

function closeDetail() {
  document.getElementById("detail").classList.add("closed");
  S.activeProt=null; S.activePd=null; S.detailMode=null; S.currentP=null;
  if (S.nglStage){S.nglStage.dispose();S.nglStage=null;}
  S.spinning=false; S.surfaceRep=null; S.hoverRep=null;
  document.querySelectorAll(".prot-row,.prot-row-ecod").forEach(el=>el.classList.remove("selected"));
  document.getElementById("sec-arch").style.display="none";
}

// ── NGL ───────────────────────────────────────────────────────────────────────
function initNGL(p, mode, pdId, container) {
  const el = container || document.getElementById("ngl-viewport");
  if (!container) {
    if (S.nglStage){S.nglStage.dispose();S.nglStage=null;}
    S.spinning=false; S.surfaceRep=null; S.hoverRep=null;
    document.getElementById("btn-spin").classList.remove("active");
    document.getElementById("btn-surf").classList.remove("active");
  }

  if (!p.pdb){el.innerHTML='<div class="ngl-no-struct">Structure not available</div>';return;}
  el.innerHTML='<div class="ngl-loading">Loading…</div>';

  try {
    const stage = new NGL.Stage(el, {backgroundColor:"white"});
    if (!container) S.nglStage=stage; else S.fsStage=stage;

    stage.loadFile(p.pdb,{ext:"pdb"}).then(comp=>{
      const loading=el.querySelector(".ngl-loading");
      if(loading)loading.remove();
      if(!container) S.nglComp=comp;

      let scheme;
      if (mode==="ecod") {
        const sel=p.tiles.map(t=>[t.c,`${t.s}-${t.e}`]); sel.push(["#aaa","*"]);
        scheme=NGL.ColormakerRegistry.addSelectionScheme(sel,"ecod_"+Date.now());
      } else if (mode==="pd"&&pdId) {
        const map=(pdByQuery.get(pdId)||[]).find(r=>r.target===p.id);
        if (map){
          const col=PD_COLOR[pdId]||"#f60";
          scheme=NGL.ColormakerRegistry.addSelectionScheme([[col,`${map.tstart}-${map.tend}`],["#ccc","*"]],"pd_"+Date.now());
        } else scheme="bfactor";
      } else scheme="bfactor";

      comp.addRepresentation("cartoon",{color:scheme,smoothSheet:true});
      comp.autoView();

      if (!container) drawSeqStrip(p, mode, pdId, stage, comp);
    }).catch(err=>{
      el.innerHTML=`<div class="ngl-no-struct">Error loading structure</div>`;
      console.error(err);
    });
  } catch(e){
    el.innerHTML=`<div class="ngl-no-struct">NGL not available – serve via HTTP</div>`;
  }
}

// ── Sequence strip ────────────────────────────────────────────────────────────
function drawSeqStrip(p, mode, pdId, stage, comp) {
  const canvas = document.getElementById("seq-canvas");
  if (!canvas) return;
  const W = canvas.parentElement.clientWidth || 400;
  canvas.width  = W;
  canvas.height = 36;
  const ctx = canvas.getContext("2d");
  const len = p.tlen || p.plen || 1;

  // Build segment list for the strip
  let segments = [];
  if (mode==="pd"&&pdId) {
    const map=(pdByQuery.get(pdId)||[]).find(r=>r.target===p.id);
    if (map) {
      if (map.tstart>1)    segments.push({s:1,e:map.tstart-1,c:"#e0e6ec",l:"Uncovered"});
      segments.push({s:map.tstart,e:map.tend,c:PD_COLOR[pdId]||"#f60",l:pdId});
      if (map.tend<len)    segments.push({s:map.tend+1,e:len,c:"#e0e6ec",l:"Uncovered"});
    } else {
      segments=[{s:1,e:len,c:"#e0e6ec",l:"No PD mapping"}];
    }
  } else {
    segments = p.tiles.map(t=>({s:t.s,e:t.e,c:t.c,l:t.l}));
  }

  function drawStrip(hoverPos) {
    ctx.clearRect(0,0,W,36);
    // background
    ctx.fillStyle="#f0f3f7"; ctx.roundRect(0,8,W,20,3); ctx.fill();
    // segments
    segments.forEach(seg => {
      const x = Math.floor((seg.s-1)/len*W);
      const w = Math.max(1, Math.ceil((seg.e-seg.s+1)/len*W));
      ctx.fillStyle=seg.c; ctx.fillRect(x,8,w,20);
    });
    // ruler ticks
    ctx.fillStyle="#bbb"; ctx.font="9px monospace";
    ctx.textAlign="left";  ctx.fillText("1",2,7);
    ctx.textAlign="right"; ctx.fillText(len,W-2,7);
    [0.25,0.5,0.75].forEach(f=>{
      const x=Math.floor(f*W);
      ctx.fillStyle="#ddd"; ctx.fillRect(x,28,1,8);
      ctx.fillStyle="#bbb"; ctx.textAlign="center"; ctx.fillText(Math.round(f*len),x,36);
    });
    // cursor
    if (hoverPos!==null) {
      const x=Math.floor((hoverPos-1)/len*W);
      ctx.strokeStyle="rgba(0,0,0,0.5)"; ctx.lineWidth=1.5;
      ctx.beginPath(); ctx.moveTo(x,4); ctx.lineTo(x,28); ctx.stroke();
    }
  }

  drawStrip(null);

  // hover
  let hoverTimer=null;
  canvas.onmousemove = e => {
    const rect=canvas.getBoundingClientRect();
    const x=(e.clientX-rect.left)*(W/rect.width);
    const pos=Math.min(len,Math.max(1,Math.round(x/W*len)));
    const seg=segments.find(s=>pos>=s.s&&pos<=s.e);
    const lbl=seg?seg.l:"—";
    drawStrip(pos);
    document.getElementById("seq-pos-label").textContent=`Position ${pos} · ${lbl}`;

    // debounced NGL highlight
    clearTimeout(hoverTimer);
    hoverTimer=setTimeout(()=>{
      if(S.nglComp===comp){
        if(S.hoverRep){comp.removeRepresentation(S.hoverRep);S.hoverRep=null;}
        try{
          S.hoverRep=comp.addRepresentation("spacefill",{sele:String(pos),color:"#FFD700",radius:1.2});
        }catch(e){}
      }
    },120);
  };
  canvas.onmouseleave = () => {
    drawStrip(null);
    document.getElementById("seq-pos-label").textContent="";
    if(S.hoverRep&&S.nglComp===comp){comp.removeRepresentation(S.hoverRep);S.hoverRep=null;}
  };
}

// ── Fullscreen ────────────────────────────────────────────────────────────────
function setupFullscreen() {
  document.getElementById("btn-fs").onclick = () => {
    const p=S.currentP;
    if(!p||!p.pdb) return;
    const fs=document.getElementById("ngl-fullscreen");
    fs.classList.add("active");
    document.getElementById("ngl-fs-title").textContent=p.id;
    if(S.fsStage){S.fsStage.dispose();S.fsStage=null;}
    const vp=document.getElementById("ngl-fs-viewport");
    vp.innerHTML="";
    initNGL(p,S.detailMode,S.activePd,vp);
  };
  document.getElementById("btn-fs-close").onclick = () => {
    document.getElementById("ngl-fullscreen").classList.remove("active");
    if(S.fsStage){S.fsStage.dispose();S.fsStage=null;}
  };
}

// ── Detail buttons ────────────────────────────────────────────────────────────
function setupDetailButtons() {
  document.getElementById("close-detail").onclick=closeDetail;
  document.getElementById("btn-spin").onclick=function(){
    S.spinning=!S.spinning; this.classList.toggle("active",S.spinning);
    if(S.nglStage)S.nglStage.setSpin(S.spinning?[0,1,0]:null);
  };
  document.getElementById("btn-reset").onclick=()=>{if(S.nglComp)S.nglComp.autoView();};
  document.getElementById("btn-surf").onclick=function(){
    if(!S.nglComp)return;
    if(S.surfaceRep){S.nglComp.removeRepresentation(S.surfaceRep);S.surfaceRep=null;this.classList.remove("active");}
    else{S.surfaceRep=S.nglComp.addRepresentation("surface",{opacity:0.3,color:"white",side:"front"});this.classList.add("active");}
  };
  document.querySelectorAll(".tab-btn").forEach(btn=>{
    btn.onclick=function(){
      document.querySelectorAll(".tab-btn").forEach(b=>b.classList.remove("active"));
      document.querySelectorAll(".tab-content").forEach(c=>c.classList.remove("active"));
      this.classList.add("active");
      document.getElementById("tab-"+this.dataset.tab).classList.add("active");
    };
  });
}

// ── Resize ────────────────────────────────────────────────────────────────────
function setupResize() {
  const h=document.getElementById("det-resize"),p=document.getElementById("detail");
  let drag=false,x0=0,w0=0;
  h.addEventListener("mousedown",e=>{drag=true;x0=e.clientX;w0=p.offsetWidth;document.body.style.cursor="ew-resize";document.body.style.userSelect="none";e.preventDefault();});
  document.addEventListener("mousemove",e=>{if(!drag)return;p.style.width=Math.min(740,Math.max(280,w0+(x0-e.clientX)))+"px";if(S.nglStage)S.nglStage.handleResize();});
  document.addEventListener("mouseup",()=>{if(drag){drag=false;document.body.style.cursor="";document.body.style.userSelect="";}});
}

// ── Tooltip ───────────────────────────────────────────────────────────────────
function setupTooltip() {
  const tip=document.getElementById("tooltip");
  document.addEventListener("mouseover",e=>{const el=e.target.closest("[data-tip]");if(el){tip.textContent=el.dataset.tip;tip.style.display="block";}});
  document.addEventListener("mousemove",e=>{tip.style.left=(e.clientX+14)+"px";tip.style.top=(e.clientY+12)+"px";});
  document.addEventListener("mouseout",e=>{if(e.target.closest("[data-tip]"))tip.style.display="none";});
}

// ── Search ────────────────────────────────────────────────────────────────────
function setupSearch() {
  const inp = document.getElementById("search-input");
  const res = document.getElementById("search-results");

  inp.addEventListener("input", () => {
    const q = inp.value.trim().toLowerCase();
    if (q.length < 2) { res.style.display = "none"; return; }

    const sections = [];

    // 1. RBP classes
    const classHits = D.classes.filter(c => c.name.toLowerCase().includes(q));
    if (classHits.length) sections.push({
      label: "RBP Classes",
      items: classHits.slice(0,5).map(c => ({
        title: "&#128193; " + esc(c.name), sub: c.n + " proteins",
        action: `jumpToClass('${esc(c.name)}')`
      }))
    });

    // 2. RBP clusters
    const clHits = [];
    D.classes.forEach(c => c.clusters.forEach(cl => {
      if (cl.id.toLowerCase().includes(q)) clHits.push({cl, cls: c.name});
    }));
    if (clHits.length) sections.push({
      label: "RBP Clusters",
      items: clHits.slice(0,5).map(({cl, cls}) => ({
        title: esc(cl.id), sub: esc(cls) + " · " + cl.proteins.length + " proteins",
        action: `jumpToCluster('${esc(cl.id)}')`
      }))
    });

    // 3. Pseudo-domain clusters
    const pdHits = PD_IDS.filter(id => id.toLowerCase().includes(q));
    if (pdHits.length) sections.push({
      label: "Pseudo-domain Clusters",
      items: pdHits.slice(0,5).map(id => ({
        title: esc(id), sub: (pdByQuery.get(id)||[]).length + " mappings",
        action: `jumpToPd('${esc(id)}')`
      }))
    });

    // 4. ECOD domains
    const ecodHits = new Map();
    D.classes.forEach(c => c.clusters.forEach(cl => cl.proteins.forEach(p =>
      p.tiles.forEach(t => {
        if (t.dk === "undetected" || t.dk === "multiple domains") return;
        if ((t.l||t.dk).toLowerCase().includes(q) && !ecodHits.has(t.dk))
          ecodHits.set(t.dk, t);
      })
    )));
    if (ecodHits.size) sections.push({
      label: "ECOD Domains",
      items: [...ecodHits.values()].slice(0,4).map(t => ({
        title: `<span style="display:inline-block;width:9px;height:9px;border-radius:2px;background:${t.c};flex-shrink:0"></span>${esc(t.l)}`,
        sub: t.dk,
        action: `jumpToEcod('${esc(t.dk)}')`
      }))
    });

    // 5. Proteins — search across many fields
    const PROT_FIELDS = [
      [p => p.id,                  "ID"],
      [p => p.gm.Phage,            "Phage"],
      [p => p.gm.Genome_name,      "Genome"],
      [p => p.gm.Genome_ID,        "Genome ID"],
      [p => p.gm.Genus,            "Genus"],
      [p => p.gm.Capsule,          "Capsule"],
      [p => p.gm.Phage_cluster,    "Phage cluster"],
      [p => p.gm.Observation,      "Observation"],
      [p => p.morph,               "Morphology"],
      [p => p.hr,                  "Host range"],
      [p => p.clid,                "Cluster"],
      [p => p.pm.PFAM_function,    "PFAM"],
      [p => p.pm.ECOD_function,    "ECOD"],
      [p => p.pm["TF-class"],      "TF class"],
    ];
    const protHits = [];
    outer: for (const [, p] of proteinMap) {
      for (const [fn, label] of PROT_FIELDS) {
        const val = fn(p);
        if (val && String(val).toLowerCase().includes(q)) {
          protHits.push({ p, label });
          if (protHits.length >= 40) break outer;
          break;
        }
      }
    }
    if (protHits.length) sections.push({
      label: "Proteins",
      items: protHits.slice(0,20).map(({p, label}) => ({
        title: `<span style="font-family:monospace">${esc(p.id)}</span>`,
        sub: esc(classOfProtein(p.id)) + " · " + label,
        action: `jumpToProtein('${esc(p.id)}')`
      }))
    });

    if (!sections.length) { res.style.display = "none"; return; }

    let html = "";
    sections.forEach(sec => {
      html += `<div class="sr-section">${sec.label}</div>`;
      sec.items.forEach(item => {
        html += `<div class="sr-item" onclick="${item.action};clearSearch()">
          <div class="sr-pid">${item.title}</div>
          <div class="sr-cls">${item.sub}</div>
        </div>`;
      });
    });
    res.innerHTML = html;
    res.style.display = "block";
  });

  document.addEventListener("click", e => {
    if (!e.target.closest(".top-search")) res.style.display = "none";
  });
}

function clearSearch(){document.getElementById("search-input").value="";document.getElementById("search-results").style.display="none";}
function classOfProtein(id){for(const c of D.classes)for(const cl of c.clusters)if(cl.proteins.some(p=>p.id===id))return c.name;return"";}
function jumpToClass(n){showView("classes");setTimeout(()=>{toggleGrp(`cls-${n}`);document.getElementById(`gb-cls-${n}`)?.scrollIntoView({behavior:"smooth",block:"start"});},80);}
function jumpToProtein(id){showView("rbp");setTimeout(()=>{openDetail(id,"basic");document.getElementById(`pr-${id}`)?.scrollIntoView({behavior:"smooth",block:"center"});},120);}
function jumpToCluster(id){showView("clusters");setTimeout(()=>{toggleGrp(`cl-${id}`);document.getElementById(`gb-cl-${id}`)?.scrollIntoView({behavior:"smooth",block:"start"});},80);}
function jumpToPd(id){showView("pd");setTimeout(()=>{toggleGrp(`pdgrp-${id}`);document.getElementById(`gb-pdgrp-${id}`)?.scrollIntoView({behavior:"smooth",block:"start"});},80);}
function jumpToEcod(dk){showView("ecod");setTimeout(()=>{toggleGrp(`ecod-${dk}`);document.getElementById(`gb-ecod-${dk}`)?.scrollIntoView({behavior:"smooth",block:"start"});},80);}

// ── Modularity network ────────────────────────────────────────────────────────
function buildModularityView() {
  // Legend
  const legEl = document.getElementById("net-legend");
  legEl.innerHTML = Object.entries(CAT_COLORS).map(([cat, col]) =>
    `<span class="cat-chip"><span class="cat-swatch" style="background:${col}"></span>${esc(cat)}</span>`
  ).join("");

  const container = document.getElementById("cy-container");
  container.innerHTML = "";
  document.getElementById("net-table-panel").style.display = "none";

  if (typeof cytoscape === "undefined") {
    container.innerHTML = '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:#aaa;font-size:13px">Cytoscape.js not loaded – serve via HTTP</div>';
    return;
  }
  if (S.cy) { S.cy.destroy(); S.cy = null; }

  const maxN    = Math.max(1, ...SM.nodes.map(n => n.n));
  const maxEdge = Math.max(1, ...SM.edges.map(e => e.n));
  const elements = [];

  SM.nodes.forEach(n => {
    const size = 22 + Math.sqrt(n.n / maxN) * 44;
    elements.push({ data: { id: n.id, label: n.id, n: n.n, size } });
  });

  SM.edges.forEach((e, i) => {
    const width = 1.5 + (e.n / maxEdge) * 6;
    elements.push({ data: {
      id: `e${i}`, source: e.s, target: e.t,
      category: e.cat, count: e.n,
      color: CAT_COLORS[e.cat] || "#aaa", width
    }});
  });

  const cy = cytoscape({
    container,
    elements,
    style: [
      { selector: "node", style: {
          "background-color": "#4a90d9",
          "label": "data(label)",
          "width": "data(size)", "height": "data(size)",
          "font-size": 9, "color": "#1a2b3c",
          "text-valign": "bottom", "text-halign": "center",
          "text-margin-y": 4, "text-wrap": "wrap", "text-max-width": 90,
          "border-width": 2, "border-color": "#2d5a8e",
          "min-zoomed-font-size": 6
      }},
      { selector: "edge", style: {
          "line-color": "data(color)", "width": "data(width)",
          "opacity": 0.55, "curve-style": "bezier",
          "loop-sweep": 100, "loop-direction": 45
      }},
      { selector: "node.highlighted", style: {
          "background-color": "#ffc107", "border-color": "#e65100",
          "border-width": 3, "z-index": 10
      }},
      { selector: "node.dimmed",  style: { "opacity": 0.15 }},
      { selector: "edge.highlighted", style: { "opacity": 0.9, "z-index": 10 }},
      { selector: "edge.dimmed",  style: { "opacity": 0.04 }}
    ],
    layout: {
      name: "cose", animate: false, randomize: true,
      nodeRepulsion: 8000, idealEdgeLength: 110, gravity: 0.6
    }
  });

  cy.on("tap", "node", function(evt) {
    const node = evt.target;
    const nid  = node.id();
    cy.elements().removeClass("highlighted dimmed");
    const hood = node.closedNeighborhood();
    hood.addClass("highlighted");
    cy.elements().not(hood).addClass("dimmed");
    S.netFilter = { type: "node", id: nid };
    showNetworkTable(nid);
  });

  cy.on("tap", "edge", function(evt) {
    const edge = evt.target;
    const src  = edge.data("source");
    const tgt  = edge.data("target");
    const cat  = edge.data("category");
    cy.elements().removeClass("highlighted dimmed");
    const sel = edge.union(edge.source()).union(edge.target());
    sel.addClass("highlighted");
    cy.elements().not(sel).addClass("dimmed");
    S.netFilter = { type: "edge", src, tgt, cat };
    showEdgeTable(src, tgt, cat);
  });

  cy.on("tap", function(evt) {
    if (evt.target === cy) {
      cy.elements().removeClass("highlighted dimmed");
      S.netFilter = null;
      document.getElementById("net-table-panel").style.display = "none";
    }
  });

  S.cy = cy;
}

// shared table renderer — shows protein-level rows
function renderSmTable(rows, title) {
  const hdr = `<div style="display:flex;align-items:center;gap:10px;margin-bottom:10px;flex-wrap:wrap">
    <strong style="font-size:12px;color:#1a2b3c">${esc(title)}</strong>
    <span style="font-size:11px;color:#888">${rows.length} protein pair${rows.length!==1?"s":""}</span>
    <button class="hdr-btn" onclick="exportNetRowsTSV()" style="margin-left:auto">&#8659; Export TSV</button>
    <button class="hdr-btn" onclick="resetNet()">&#10005;</button>
  </div>`;

  const cols = `<thead><tr>
    <th>Query</th><th>Target</th><th>Category</th>
    <th title="% identity">%id</th><th title="Alignment length">Aln</th>
    <th>Q region</th><th>T region</th>
    <th title="Query coverage">qcov</th><th title="Target coverage">tcov</th>
    <th>E-value</th>
  </tr></thead>`;

  let body = "<tbody>";
  rows.forEach(r => {
    const dot = `<span style="width:8px;height:8px;border-radius:50%;flex-shrink:0;display:inline-block;background:${CAT_COLORS[r.cat]||"#aaa"}"></span>`;
    const ev  = r.ev === 0 ? "0" : r.ev.toExponential(1);
    body += `<tr>
      <td style="font-family:monospace;font-size:10px">${esc(r.q)}</td>
      <td style="font-family:monospace;font-size:10px">${esc(r.t)}</td>
      <td><span style="display:inline-flex;align-items:center;gap:4px">${dot}${esc(r.cat)}</span></td>
      <td>${r.pid}</td><td>${r.al}</td>
      <td style="white-space:nowrap">${r.qs}–${r.qe}</td>
      <td style="white-space:nowrap">${r.ts}–${r.te}</td>
      <td>${(r.qcov*100).toFixed(0)}%</td>
      <td>${(r.tcov*100).toFixed(0)}%</td>
      <td style="font-family:monospace;font-size:10px">${ev}</td>
    </tr>`;
  });
  body += "</tbody>";

  const panel = document.getElementById("net-table-panel");
  panel.innerHTML = hdr + `<table class="net-table">${cols}${body}</table>`;
  panel.style.display = "block";
}

function showNetworkTable(nodeId) {
  const rows = SM_ROWS.filter(r => r.qc === nodeId || r.tc === nodeId);
  rows.sort((a, b) => a.cat.localeCompare(b.cat) || b.pid - a.pid);
  renderSmTable(rows, nodeId);
}

function showEdgeTable(src, tgt, cat) {
  const rows = SM_ROWS.filter(r =>
    r.cat === cat &&
    ((r.qc === src && r.tc === tgt) || (r.qc === tgt && r.tc === src))
  );
  rows.sort((a, b) => b.pid - a.pid);
  const isSelf  = src === tgt;
  const partner = isSelf ? `${esc(src)} (self)` : `${esc(src)} ↔ ${esc(tgt)}`;
  renderSmTable(rows, `${partner} · ${esc(cat)}`);
}

function resetNet() {
  if (S.cy) S.cy.elements().removeClass("highlighted dimmed");
  S.netFilter = null;
  document.getElementById("net-table-panel").style.display = "none";
}

function exportNetRowsTSV() {
  const f = S.netFilter;
  let rows;
  if (!f) {
    rows = SM_ROWS;
  } else if (f.type === "node") {
    rows = SM_ROWS.filter(r => r.qc === f.id || r.tc === f.id);
  } else {
    rows = SM_ROWS.filter(r =>
      r.cat === f.cat &&
      ((r.qc === f.src && r.tc === f.tgt) || (r.qc === f.tgt && r.tc === f.src))
    );
  }
  const lines = ["query\ttarget\tquery_class\ttarget_class\tcategory\tpident\taln_len\tqstart\tqend\ttstart\ttend\tqcov\ttcov\tevalue"];
  rows.forEach(r => lines.push(
    `${r.q}\t${r.t}\t${r.qc}\t${r.tc}\t${r.cat}\t${r.pid}\t${r.al}\t${r.qs}\t${r.qe}\t${r.ts}\t${r.te}\t${r.qcov}\t${r.tcov}\t${r.ev}`
  ));
  const blob = new Blob([lines.join("\n")], { type: "text/tab-separated-values" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  const tag = !f ? "all" : f.type === "node" ? f.id : `${f.src}-${f.tgt}`;
  a.download = `rbp-modularity-${tag}-${new Date().toISOString().slice(0,10)}.tsv`;
  a.click(); URL.revokeObjectURL(a.href);
}

// ── Helpers ───────────────────────────────────────────────────────────────────
function makeTable(obj,labels){
  let rows="";
  for(const[k,lbl]of Object.entries(labels)){
    const v=obj[k]||""; if(!v||v==="nan"||v==="NaN")continue;
    let vh=esc(v);
    if(k==="PMID")vh=v.split(/[,;]\s*/).map(t=>{t=t.trim();if(!t)return"";return /^\d+$/.test(t)?`<a class="pmid-link" href="https://pubmed.ncbi.nlm.nih.gov/${t}" target="_blank">${t}</a>`:esc(t);}).join(", ");
    rows+=`<tr><td>${esc(lbl)}</td><td>${vh}</td></tr>`;
  }
  return `<table class="meta-table"><tbody>${rows}</tbody></table>`;
}
function esc(s){return String(s||"").replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;").replace(/'/g,"&#39;");}
</script>
</body>
</html>
"""

# ── 11. Inject + write ─────────────────────────────────────────────────────────
html_out = (HTML
    .replace("__ATLAS_DATA__", data_json)
    .replace("__PD_DATA__",    pd_json)
    .replace("__SM_DATA__",    sm_json)
    .replace("__SM_ROWS__",    sm_rows_json)
)
OUT.write_text(html_out, encoding="utf-8")
print(f"Written: {OUT}  ({OUT.stat().st_size//1024} KB)")

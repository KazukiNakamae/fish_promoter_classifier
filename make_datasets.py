#!/usr/bin/env python3
# make_datasets.py
import argparse, os, sys, glob, re, csv, json
from pathlib import Path
from collections import defaultdict
import random

def list_fasta(paths):
    out = []
    for p in paths:
        P = Path(p)
        if P.is_dir():
            out += sorted([str(x) for x in P.glob("*.fa")] + [str(x) for x in P.glob("*.fasta")])
        elif P.is_file():
            out.append(str(P))
        else:
            sys.exit(f"[ERROR] Not found: {p}")
    if not out:
        sys.exit("[ERROR] No FASTA files found from inputs.")
    return out

def read_fasta(fa_path):
    with open(fa_path) as f:
        name=None; buf=[]
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = line[1:].split()[0]
                buf=[]
            else:
                buf.append(line)
        if name is not None:
            yield name, "".join(buf)

def load_species_map(species_map_path):
    # CSV (file_basename_without_ext,species) or JSON {"substr":"Species", ...}
    if not species_map_path: return {}
    p = Path(species_map_path)
    m = {}
    if p.suffix.lower() == ".json":
        with open(p) as f:
            m = json.load(f)
        # JSON は { "substr": "Species" } として扱い、ファイル名に substr が含まれれば species とする
        return {"__substr__": m}
    else:
        # CSV: basename_without_ext,species
        with open(p) as f:
            for row in csv.reader(f):
                if not row: continue
                m[row[0]] = row[1]
        return m

def infer_species(basename, spmap):
    if not spmap: return None
    if "__substr__" in spmap:
        for substr, spp in spmap["__substr__"].items():
            if substr in basename:
                return spp
        return None
    return spmap.get(basename)

def main():
    ap = argparse.ArgumentParser(description="Create DNABERT-2 train/dev/test CSV from positive/negative FASTA files.")
    ap.add_argument("--pos", nargs="+", required=True, help="Positive FASTA files or directories (one or many)")
    ap.add_argument("--neg", nargs="+", required=True, help="Negative FASTA files or directories (one or many)")
    ap.add_argument("--out_dir", required=True, help="Output directory for CSVs")
    ap.add_argument("--expect_len", type=int, default=700, help="Expected sequence length; others will be dropped")
    ap.add_argument("--drop_N", action="store_true", help="Drop sequences containing N (recommended)")
    ap.add_argument("--uppercase", action="store_true", help="Uppercase sequences")
    ap.add_argument("--test_size", type=float, default=0.10, help="Fraction for test split")
    ap.add_argument("--dev_size", type=float, default=0.10, help="Fraction for dev split (of total)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--stratify", choices=["none","label","species","species_label"], default="species_label",
                    help="Stratification key")
    ap.add_argument("--species_map", default=None, help="CSV (basename,species) or JSON {substr: species}")
    ap.add_argument("--keep_species", action="store_true", help="Keep species column in output CSVs (third column)")
    args = ap.parse_args()

    random.seed(args.seed)
    os.makedirs(args.out_dir, exist_ok=True)
    spmap = load_species_map(args.species_map)

    # gather fasta paths
    pos_fas = list_fasta(args.pos)
    neg_fas = list_fasta(args.neg)

    rows = []  # (seq, label, species)
    for path, label in [(p,1) for p in pos_fas] + [(n,0) for n in neg_fas]:
        bn = os.path.basename(path)
        base = os.path.splitext(bn)[0]
        species = infer_species(base, spmap)
        for name, seq in read_fasta(path):
            s = seq
            if args.uppercase: s = s.upper()
            s = re.sub(r'[^ACGTNacgtn]', '', s)
            if args.drop_N and ("N" in s.upper()):
                continue
            if args.expect_len and len(s) != args.expect_len:
                continue
            rows.append( (s.upper(), label, species) )

    if not rows:
        sys.exit("[ERROR] No sequences passed filters. Check input and --expect_len/--drop_N.")

    # stratified split
    def key(row):
        seq,label,sp = row
        if args.stratify == "none": return "all"
        if args.stratify == "label": return f"L{label}"
        if args.stratify == "species": return f"S{sp}"
        return f"S{sp}_L{label}"  # species_label

    bykey = defaultdict(list)
    for r in rows:
        bykey[key(r)].append(r)

    # compute sizes
    N = len(rows)
    n_test = int(round(N * args.test_size))
    n_dev  = int(round(N * args.dev_size))
    n_train = N - n_test - n_dev
    if n_train <= 0:
        sys.exit("[ERROR] Split sizes too small. Adjust --test_size/--dev_size.")

    # sample per stratum with proportional allocation
    def stratified_sample(pool, k):
        if k <= 0: return [], pool
        take = min(k, len(pool))
        idxs = list(range(len(pool)))
        random.shuffle(idxs)
        sel = [pool[i] for i in idxs[:take]]
        rem = [pool[i] for i in idxs[take:]]
        return sel, rem

    # flatten helpers
    train=[]; dev=[]; test=[]
    # determine per-key quotas
    quotas_test = {}
    quotas_dev = {}
    for k,v in bykey.items():
        f = len(v)/N
        quotas_test[k] = int(round(n_test * f))
        quotas_dev[k]  = int(round(n_dev  * f))

    # per-key sample
    leftovers=[]
    for k,v in bykey.items():
        random.shuffle(v)
        sel_test, rest = stratified_sample(v, quotas_test[k])
        sel_dev,  rest = stratified_sample(rest, quotas_dev[k])
        test += sel_test; dev += sel_dev; leftovers += rest
    train = leftovers

    # write CSVs (DNABERT-2 expects header: sequence,label)
    def write_csv(path, items):
        with open(path, "w", newline="") as g:
            w = csv.writer(g)
            if args.keep_species:
                w.writerow(["sequence","label","species"])
                for s,l,sp in items:
                    w.writerow([s,l,sp if sp is not None else ""])
            else:
                w.writerow(["sequence","label"])
                for s,l,sp in items:
                    w.writerow([s,l])

    write_csv(os.path.join(args.out_dir, "train.csv"), train)
    write_csv(os.path.join(args.out_dir, "dev.csv"),   dev)
    write_csv(os.path.join(args.out_dir, "test.csv"),  test)

    print(f"[DONE] Wrote CSVs to {args.out_dir}")
    print(f" train: {len(train)}; dev: {len(dev)}; test: {len(test)}")
    print(" Tip: DNABERT-2 expects three CSVs with header 'sequence,label'. See official docs.", file=sys.stderr)

if __name__ == "__main__":
    main()

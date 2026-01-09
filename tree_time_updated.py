#!/usr/bin/env python3
"""
If the proportion of non-empty date values in metadata.tsv.gz is high enough,
then get the RefSeq length from output_stats.tsv, make a newick file with
branch lengths scaled as substitutions per site, make a dates.csv file with
added -XX for missing month and/or day, and run treetime to reroot the tree.
Then apply the same rerooting to viz.pb.gz.
"""

import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
import treeswift

import viral_usher_trees
import alter_gbff

default_min_real_dates = 0.8


def get_dates(subdir_path):
    """Scan metadata.tsv.gz, return a dict of names -> dates and the proportion of real date values to total count"""
    name_to_date = {}
    with gzip.open(subdir_path + "/metadata.tsv.gz", "rt") as f:
        header = f.readline().split('\t')
        for idx, field in enumerate(header):
            header[idx] = field.strip()
        name_idx = header.index('strain')
        if name_idx != 0:
            print(f"Uh-oh, name_idx is {name_idx}", file=sys.stderr)
            sys.exit(1)
        date_idx = header.index('date')
        if date_idx < 0:
            print(subdir_path + "/metadata.tsv.gz" + " does not have date column", file=sys.stderr)
            sys.exit(1)
        line_count = 0
        real_date_count = 0
        for line in f:
            fields = line.rstrip().split("\t")
            name = fields[name_idx]
            date = fields[date_idx]
            if re.match('^[0-9]{4}', date):
                name_to_date[name] = date
                real_date_count += 1
            line_count += 1
    return name_to_date, (real_date_count / line_count)


def scale_branch_lengths(subdir_path, newick_out, refseq_len):
    """Use treeswift (included with taxoniumtools) to scale branch lengths as expected by treetime"""
    tree = treeswift.read_tree_newick(subdir_path + "/viz.nwk.gz")
    tree.scale_edges(1.0 / refseq_len)
    tree.write_tree_newick(subdir_path + "/" + newick_out)


def make_dates_csv(subdir_path, dates_out, name_to_date):
    """Write a CSV file of names and dates, fillin in missing month and/or day with -XX"""
    with open(subdir_path + "/" + dates_out, "w") as f:
        f.write(",".join(["name", "date"]) + "\n")
        for name, date in name_to_date.items():
            if re.match('^[0-9]{4}$', date):
                date += "-XX-XX"
            elif re.match('^[0-9]{4}-[0-9]{2}$', date):
                date += "-XX"
            f.write(",".join([name, date]) + "\n")


def get_refseq_len(subdir_path):
    """Get refseq_length value from output_stats.tsv"""
    with open(subdir_path + "/output_stats.tsv", "r", newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        if "ref_length" not in fieldnames:
            print("output_stats.tsv file does not have ref_length column", file=sys.stderr)
            sys.exit(1)
        for row in reader:
            ref_length = int(row["ref_length"])
    return ref_length


def run_treetime(path, min_real_dates):
    """Format input for treetime clock, run it and apply the same rooting to viz.pb.gz"""
    tree = path.strip("/") + "/optimized.pb.gz"
    path = path.strip("/")

    #subdir_path = viral_usher_trees.trees_dir + "/" + tree
    #print(subdir_path)
    #if not os.path.isdir(subdir_path):
    #    print(f"{tree} does not have a subdirectory {tree}, check spelling", file=sys.stderr)
    #    sys.exit(1)
    if not os.path.isfile(tree):
        print(f"Tree file {tree} not found, check spelling", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(path + "/metadata.tsv.gz") or not os.path.exists(path + "/output_stats.tsv"):
        print(f"Expected files not found in {path}, has the tree been built?", file=sys.stderr)
        sys.exit(1)
    
    name_to_date, real_dates_proportion = get_dates(path)
    if real_dates_proportion < min_real_dates:
        print(f"Tree {tree} has too low a proportion of dates ({real_dates_proportion:.2f} < {min_real_dates:.2f}), not running treetime. ")
        sys.exit(0)
    refseq_len = get_refseq_len(path)
    #made it here
    scale_branch_lengths(path, "viz.scaled.nwk", refseq_len)
    make_dates_csv(path, "dates.csv", name_to_date)
    command = ["treetime", "clock",
               "--sequence-length", str(refseq_len),
               "--tree",  path + "/viz.scaled.nwk",
               "--dates", path + "/dates.csv",
               "--outdir", path + "/treetime_out"]
    try:
        subprocess.run(command, check=True)
    except Exception as e:
        print(f"treetime command ({' '.join(command)}) failed: {e}", file=sys.stderr)
        sys.exit(1)


def get_oldest_node(tree):
    """Extract the oldest internal node date from treetime's output file rtt.csv"""
    # oldest_node=$(grep ^node_ treetime_out/rtt.csv | sort -t, -k2n | head -1 | cut -d, -f 1)
    tree= tree.strip("/")
    rtt_csv_path = "/".join([tree, "treetime_out", "rtt.csv"])
    min_date = None
    oldest_node = None
    with open(rtt_csv_path) as f:
        for line in f:
            columns = line.split(",")
            if columns[0].startswith("node_"):
                date = float(columns[1].strip())
                if min_date is None or date < min_date:
                    min_date = date
                    oldest_node = columns[0].strip()
    if oldest_node is None:
        print(f"Failed to get oldest node from {rtt_csv_path}", file=sys.stderr)
        sys.exit(1)
    return oldest_node


def get_refseq_acc(config_path):
    """Return the value of refseq_acc from config.toml"""
    refseq_acc = None
    with open(config_path) as f:
        for line in f:
            if line.startswith("refseq_acc = "):
                refseq_acc = line.split("'")[1]
    if refseq_acc is None:
        print(f"Failed to find refseq_acc in {config_path}", file=sys.stderr)
        sys.exit(1)
    return refseq_acc


def reroot_tree(tree, oldest_node):
    """Reroot viz.pb.gz to oldest_node using matUtils and return path to rerooted .pb.gz."""
    tree = tree.strip("/")
    input_path = "/".join([tree, "viz.pb.gz"])
    rerooted_tree_path = "/".join([tree, "timetree_rerooted.pb.gz"])
    config_path = "/".join([tree, "config.toml"])
    refseq_acc = get_refseq_acc(config_path)
    ref_path = "/".join([tree, refseq_acc + ".fa"])
    modified_ref_path = "/".join([tree, "treetime_rerooted_" + refseq_acc + ".fa"])
    command = ["matUtils", "extract",
               "-i", input_path,
               "--reroot", oldest_node,
               "--input-fasta", ref_path,
               "--write-reroot-reference", modified_ref_path,
               "-o", rerooted_tree_path]
    try:
        subprocess.run(command, check=True)
    except Exception as e:
        print(f"matUtils command ({' '.join(command)}) failed: {e}", file=sys.stderr)
        sys.exit(1)
    gbff_path = "/".join([tree, refseq_acc + ".gbff"])
    rerooted_gbff_path = "/".join([tree, "treetime_rerooted_" + refseq_acc + ".gbff"])
    alter_gbff.alter_gbff_file(gbff_path, refseq_acc, ref_path, rerooted_gbff_path)
    return rerooted_tree_path, rerooted_gbff_path


def get_columns_string(tsv_path):
    """Return a comma-separated string listing columns of TSV file."""
    # columns=$(gunzip -c metadata.tsv.gz | head -1 | sed -re 's/\t/,/g')
    with gzip.open(tsv_path, "rt") if tsv_path.endswith(".gz") else open(tsv_path) as f:
        header = f.readline().split('\t')
        columns = ",".join([col.strip() for col in header])
        return columns
    print(f"Failed to get columns from {tsv_path}", file=sys.stderr)
    sys.exit(1)


def make_taxonium(tree, rerooted_tree_path, rerooted_gbff_path):
    """Run usher_to_taxonium on rerooted_tree"""
    tree = tree.strip("/")
    metadata_path = "/".join([tree, "metadata.tsv.gz"])
    columns = get_columns_string(metadata_path)
    taxonium_jsonl_path = "/".join([tree, "timetree_rerooted.jsonl.gz"])
    command = ["usher_to_taxonium",
               "-i", rerooted_tree_path,
               "-m", metadata_path,
               "--genbank", rerooted_gbff_path,
               "-c", columns,
               "--title", "Treetime-rerooted " + tree,
               "-o", taxonium_jsonl_path]
    try:
        subprocess.run(command, check=True)
    except Exception as e:
        print(f"usher_to_taxonium command ({' '.join(command)}) failed: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Use tree metadata TSV and trees/*/output_stats.tsv to make summary TSV")
    parser.add_argument('-t', '--virus_path', required=True,
                        help="Path to virus dir")
    parser.add_argument('-m', '--min_real_dates',
                        help=f"Minimum proportion of dates in metadata.tsv.gz that have real values (default: {default_min_real_dates})")
    args = parser.parse_args()
    #viral_usher_trees.check_top_level_dir()
    min_real_dates = args.min_real_dates if args.min_real_dates else default_min_real_dates
    #tree = args.virus_path.strip("/") + "/optimized.pb.gz"
    #print(tree)
    run_treetime(args.virus_path, min_real_dates)
    oldest_node = get_oldest_node(args.virus_path)
    rerooted_tree, rerooted_gbff = reroot_tree(args.virus_path, oldest_node)
    make_taxonium(args.virus_path, rerooted_tree, rerooted_gbff)


if __name__ == "__main__":
    main()
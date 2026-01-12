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

#import viral_usher_trees

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
            print(name)
            print(date)
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
            print(name)
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
        if "refseq_length" not in fieldnames:
            return -1
        for row in reader:
            refseq_length = int(row["refseq_length"])
    return refseq_length


def run_treetime(virus, directory, min_real_dates):
    """Format input for treetime clock, run it and apply the same rooting to viz.pb.gz"""
    subdir_path = directory + "/" + virus
    if not os.path.isdir(directory):
        print(f"{directory} does not have a subdirectory {virus}, check spelling", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(subdir_path + "/metadata.tsv.gz") or not os.path.exists(subdir_path + "/output_stats.tsv"):
        print(f"Expected files not found in {subdir_path}, has the tree been built?", file=sys.stderr)
        sys.exit(1)
    name_to_date, real_dates_proportion = get_dates(subdir_path)
    if real_dates_proportion < min_real_dates:
        with open(subdir_path + "/treetime.log", "w") as outfile:
            outfile.write(f"Tree {virus} has too low a proportion of dates ({real_dates_proportion:.2f} < {min_real_dates:.2f}), not running treetime.\n")

        print(f"Tree {virus} has too low a proportion of dates ({real_dates_proportion:.2f} < {min_real_dates:.2f}), not running treetime. ")
        sys.exit(0)
    refseq_len = get_refseq_len(subdir_path)
    scale_branch_lengths(subdir_path, "viz.scaled.nwk", refseq_len)
    make_dates_csv(subdir_path, "dates.csv", name_to_date)
    command = ["treetime", "clock",
               "--sequence-length", str(refseq_len),
               "--tree",  subdir_path + "/viz.scaled.nwk",
               "--dates", subdir_path + "/dates.csv",
               "--outdir", subdir_path + "/treetime_out"]
    '''
    with open(subdir_path + "/treetime.log", "w") as outfile:
        try:
            subprocess.run(command, check=True, stdout=outfile, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print(f"treetime command ({' '.join(command)}) failed with return code {e.returncode}", file=outfile)
            sys.exit(1)
        finally:
            outfile.write("\n[Snakemake wrapper] treetime finished (check above for success/failure)\n")
    '''
    '''
    with open(subdir_path + "/treetime.log", "w") as outfile:
        try:
            subprocess.run(command, check=True, stdout=outfile)
        except Exception as e:
            print(f"treetime command ({" ".join(command)}) failed: {e}", file=outfile)
            sys.exit(1)
    '''
    with open(subdir_path + "/treetime.log", "w") as outfile:
        try:
            subprocess.run(command, check=True, stdout=outfile, stderr=subprocess.STDOUT)
            outfile.write("\n[OK] treetime finished successfully\n")
        except subprocess.CalledProcessError as e:
            outfile.write(f"\n[ERROR] treetime command ({' '.join(command)}) failed with return code {e.returncode}\n")
            # DO NOT sys.exit(1) — just log the failure
        
        #print(" TODO IMPLEMENT ME figure out how to apply the same rerooting to viz.pb.gz, i.e. which node did treetime choose?")


def main():
    parser = argparse.ArgumentParser(description="Use tree metadata TSV and trees/*/output_stats.tsv to make summary TSV")
    parser.add_argument('-v', '--virus', required=True,
                        help="Tree name (must match a name from the tree_name column of tree_metadata.tsv)")
    parser.add_argument('-m', '--min_real_dates',
                        help=f"Minimum proportion of dates in metadata.tsv.gz that have real values (default: {default_min_real_dates})")
    parser.add_argument('-d', '--directory', help="Directory containing trees", required=True)
    args = parser.parse_args()
    #viral_usher_trees.check_top_level_dir()
    min_real_dates = args.min_real_dates if args.min_real_dates else default_min_real_dates
    run_treetime(args.virus, args.directory, min_real_dates)


if __name__ == "__main__":
    main()

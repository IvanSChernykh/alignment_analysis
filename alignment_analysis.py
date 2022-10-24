import sys as ss
import os
import pandas as pd
import copy as cp

def alignment_split(indir, outdir):
    try:
        os.mkdir(f"{outdir}\\alan_result")
    except FileExistsError:
        pass
    try:
        os.mkdir(f"{outdir}\\alan_result\\align_split")
    except FileExistsError:
        pass
    with open(f"{indir}\\pangenome\\pangenome.bs") as pangenome:
        blocks_list = []
        for line in pangenome:
            if line.startswith(">"):
                line_split = line.split()
                block = line_split[1][6:]
                if block not in blocks_list:
                    blocks_list.append(block)
                    with open(f"{outdir}\\alan_result\\align_split\\{block}.fasta", "w") as block_file:
                        block_file.write(line)
                else:
                    with open(f"{outdir}\\alan_result\\align_split\\{block}.fasta", "a") as block_file:
                        block_file.write(f"\n{line}")
            else:
                with open(f"{outdir}\\alan_result\\align_split\\{block}.fasta", "a") as block_file:
                    block_file.write(line.strip())

def alignment_analysis(indir, outdir, over, split_dir):
    
    def group_finder(column, group):
        result = []
        for elem in column:
            result.append(group in elem.split())
        return(result)
    
    table = pd.read_csv(f"{indir}\\genes\\partition-ungrouped.tsv", sep = "\t")
    table = table.iloc[:, :-1]
    table["end"] = [False]*len(table)
    table_copy = table
    for gene, loc_table in table.groupby("gene"):
        loc_table = loc_table.sort_values("gene_block_stop", ascending = False)
        table_copy.loc[loc_table.iloc[0].name, "end"] = True
    table = table_copy
    table["group"] = [""]*len(table)
    group_id = 0
    table_copy = table
    for block, loc_table in table[table["gene_block_start"] == 0].groupby("npg_block"):
        loc_table = loc_table.sort_values("npg_block_min")
        indexes = []
        for index in loc_table.index:
            if loc_table.loc[index, "npg_block_max"] - loc_table.loc[index, "npg_block_min"] < over:
                table_copy.loc[index, "group"] += f" {group_id}"
                group_id += 1
            else:
                indexes.append(index)
                potential = loc_table.loc[indexes]
                end = min(potential["npg_block_max"])
                success = False
                while (len(indexes) > 1) and (not success):
                    if end - loc_table.loc[index, "npg_block_min"] >= over:
                        success = True
                    else:
                        table_copy.loc[indexes[:-1], "group"] += f" {group_id}"
                        group_id += 1
                        new_indexes = []
                        for i in range(len(indexes) - 1):
                            if potential.loc[indexes[i], "npg_block_max"] != end:
                                new_indexes.append(indexes[i])
                        new_indexes.append(indexes[-1])
                        indexes = cp.copy(new_indexes)
                        potential = loc_table.loc[indexes]
                        end = min(potential["npg_block_max"])
        if len(indexes) != 0:
            table_copy.loc[indexes, "group"] += f" {group_id}"
            group_id += 1
    table = table_copy
    for gene, loc_table in table.groupby("gene"):
        loc_table = loc_table.sort_values("gene_block_start")
        if len(loc_table[loc_table["gene_block_start"] == 0]) == 1:
            indexes = loc_table.index
            for i in range(1, len(indexes)):
                table_copy.loc[indexes[i], "group"] = loc_table.loc[indexes[0], "group"]
    table = table_copy
    table["New block start"] = [-1]*len(table)
    group_table = {"Group": [], "Mistake start": [], "Ori warning": [], "Len warning": []}
    for i in range(group_id):
        group_table["Group"].append(i)
        mistake_start = "No"
        ori_warning = False
        len_warning = False
        loc_table = table[group_finder(table["group"], str(i))]
        loc_table_start = loc_table[loc_table["gene_block_start"] == 0]
        block = loc_table_start.iloc[0, 4]
        loc_table_start = loc_table_start.sort_values("npg_block_min")
        group_start_first = loc_table_start.iloc[0, 5]
        seq_first = loc_table_start.iloc[0, 0]
        gene_first = loc_table_start.iloc[0, 3]
        ori_first = loc_table_start.iloc[0, 7]
        if loc_table_start.iloc[0, 1] <= loc_table_start.iloc[0, 2]:
            seq_start_first = loc_table_start.iloc[0, 1]
            seq_stop_first = loc_table_start.iloc[0, 2]
            ori_first = True
        else:
            seq_start_first = loc_table_start.iloc[0, 2]
            seq_stop_first = loc_table_start.iloc[0, 1]
            ori_first = False
        indexes = loc_table_start.index[1:]
        for index in indexes:
            if loc_table_start.loc[index, "npg_block_ori"] != ori_first:
                ori_warning = True
            if group_start_first != loc_table_start.loc[index, "npg_block_min"]:
                if (loc_table_start.loc[index, "npg_block_min"] - group_start_first)%3 == 0:
                    len_warning = True
                sequence_first = "first"
                sequence_second = "second"
                seq_second = loc_table_start.loc[index, "sequence"]
                group_start_second = loc_table_start.loc[index, "npg_block_min"]
                if loc_table_start.loc[index, "sequence_start"] <= loc_table_start.loc[index, "sequence_stop"]:
                    seq_start_second = loc_table_start.loc[index, "sequence_start"]
                    seq_stop_second = loc_table_start.loc[index, "sequence_stop"]
                    ori_second = True
                else:
                    seq_start_second = loc_table_start.loc[index, "sequence_stop"]
                    seq_stop_second = loc_table_start.loc[index, "sequence_start"]
                    ori_second = False
                with open(f"{split_dir}\\{block}.fasta") as block_file:
                    succ_first = False
                    succ_second = False
                    line = block_file.readline().strip()
                    while (line != "") and (not (succ_first and succ_second)):
                        if line.startswith(">"):
                            line = line[1:].split()[0].split("_")
                            seq = line[0]
                            if int(line[1]) <= int(line[2]):
                                seq_start = int(line[1])
                                seq_stop = int(line[2])
                            else:
                                seq_start = int(line[2])
                                seq_stop = int(line[1])
                            if (seq == seq_first) and (seq_start <= seq_start_first) and (seq_stop >= seq_stop_first):
                                succ_first = True
                                sequence_first = block_file.readline().strip()
                                if ori_first:
                                    sequence_first = sequence_first[seq_start_first - seq_start:group_start_second - group_start_first + seq_start_first - seq_start + 1]
                                else:
                                    sequence_first = sequence_first[seq_start - seq_start_first:group_start_second - group_start_first + seq_start - seq_start_first + 1]
                            if (seq == seq_second) and (seq_start <= seq_start_second) and (seq_stop >= seq_stop_second):
                                succ_second = True
                                sequence_second = block_file.readline().strip()
                                if ori_second:
                                    sequence_second = sequence_second[seq_start_second - seq_start - group_start_second + group_start_first:seq_start_second - seq_start + 1]
                                else:
                                    sequence_second = sequence_second[seq_start - seq_start_second - group_start_second + group_start_first:seq_start - seq_start_second + 1]
                        line = block_file.readline().strip()
                if sequence_first == sequence_second:
                    mistake_start = "Yes"
                    table.loc[index, "New block start"] = group_start_first
        group_table["Ori warning"].append(ori_warning)
        group_table["Mistake start"].append(mistake_start)
        group_table["Len warning"].append(len_warning)
    group_table = pd.DataFrame(group_table)
    table.to_excel(f"{outdir}\\alan_result\\partition_remastered.xlsx")
    group_table.to_excel(f"{outdir}\\alan_result\\group_result.xlsx")
    start_number = len(group_table[group_table["Mistake start"] == "Yes"])
    ori_number = len(group_table[group_table["Ori warning"]])
    len_number = len(group_table[group_table["Len warning"]])
    remastered = len(table[table["New block start"] != -1])
    n_genes = len(table)
    with open(f"{outdir}\\alan_result\\general_result.txt", "w") as general:
        general.write(f"Total number of groups: {group_id}\n")
        general.write(f"Total number of groups with mistakes in start: {start_number}\n")
        general.write(f"Total number of groups with warning in gene orientation: {ori_number}\n")
        general.write(f"Total number of groups with warning in gene length: {len_number}\n")
        general.write(f"Total number of genes: {n_genes}\n")
        general.write(f"Total number of reannotated genes: {remastered}\n")
    
if "-h" in ss.argv:
    print("Alignment analysis, version 1.1.0")
    print("Description of parameters (all are optional):")
    print("-h - type this message and exit")
    print("-i - input directory (default - current)")
    print("-o - output directory (default - current)")
    print("-l - overlap in nucleotides between genes to form group (default - 90)")
    print("-a - if given, the alignment split would not be performed")
else:
    try:
        indir_index = ss.argv.index("-i")
    except ValueError:
        indir = os.getcwd()
    else:
        indir = ss.argv[indir_index + 1]
    try:
        outdir_index = ss.argv.index("-o")
    except ValueError:
        outdir = os.getcwd()
    else:
        outdir = ss.argv[outdir_index + 1]
        try:
            os.mkdir(outdir)
        except FileExistsError:
            pass
    if "-a" not in ss.argv:
        alignment_split(indir, outdir)
    align = f"{outdir}\\alan_result\\align_split"
    try:
        over_index = ss.argv.index("-l")
    except ValueError:
        over = 90
    else:
        over = int(ss.argv[over_index + 1])
    alignment_analysis(indir, outdir, over, align)
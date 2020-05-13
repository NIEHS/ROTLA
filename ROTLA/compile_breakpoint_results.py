#!/usr/bin/env python

import sys

def compile_breakpoints(input_files, output_file):

    breakpoint_dict = dict()
    id_list = []
    totals_dict = dict()

    def getBreaks(file_name, file_id):
        with open(file_name) as f:
            next(f)
            for line in f:
                breakpoint_0, breakpoint_1, count = line.strip().split("\t")
                count = int(count)
                if (breakpoint_0, breakpoint_1) in breakpoint_dict:
                    breakpoint_dict[(breakpoint_0, breakpoint_1)][file_id] = count
                else:
                    breakpoint_dict[(breakpoint_0, breakpoint_1)] = {file_id:count}

    with open(input_files) as file_list:
        for line in file_list:
            file_name, file_id = line.strip().split()
            getBreaks(file_name, file_id)
            id_list.append(file_id)

    ## ADD ZEROES
    for breakpoint in breakpoint_dict:
        for file_id in id_list:
            if file_id not in breakpoint_dict[breakpoint]:
                breakpoint_dict[breakpoint][file_id] = 0

    ## GET TOTALS
    for breakpoint in breakpoint_dict:
        total = 0
        for file_id in id_list:
            total += breakpoint_dict[breakpoint][file_id]
        totals_dict[breakpoint] = total

    ## SORT BREAKPOINT LIST BY TOTALS
    sorted_breakpoint = sorted(totals_dict, key=lambda k: totals_dict[k], reverse=True)

    with open(output_file, "w") as OUTPUT:
        OUTPUT.write("\t")
        for file_id in id_list:
            OUTPUT.write("\t" + file_id)
        OUTPUT.write("\n")

        for breakpoint in sorted_breakpoint:
            OUTPUT.write(breakpoint[0] + "\t" + breakpoint[1])
            for file_id in id_list:
                OUTPUT.write("\t" + str(breakpoint_dict[breakpoint][file_id]))
            OUTPUT.write("\n")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.stdout.write("Usage: " + sys.argv[0] + "\n           <List of breakpoint files>\n           <Output file name>\n")
        exit()
    else:
        compile_breakpoints(sys.argv[1], sys.argv[2])

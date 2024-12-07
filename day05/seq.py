import sys
import os.path

def main():
    # validate the argv inputs
    argv = sys.argv

    validate_input(argv)
    
    all = []
    dict_all = {"A": 0, "T": 0, "C": 0, "G": 0, "Unknown": 0, "Total": 0}

    # go over each file given and output the results
    for path in argv[1:]:
        # extraction the seq from each file
        seq = read_file(path)

        # saving the analysis results in a list
        all.append(analyze_seq(seq))

    for i, result_dict in enumerate(all):
        # printing the file name
        print(argv[i + 1])

        # printing the analysis results
        print_result(result_dict)

        # updating the all the results to a single dict
        for nuc,count in result_dict.items():
             dict_all[nuc] += count
    
    # printing the analysis of all the files combined
    print("All")
    print_result(dict_all)


def validate_input(argv):
    # check for file paths
    if len(argv) < 2:
        exit("Usage: python seq.py filepath1.txt filepath2.txt")

    # check that all files are .txt files
    for path in argv[1:]:
        if path[-4:] != '.txt':
            exit(f"Usage: '{path}' is not a .txt file")
    
        # validate that the file exists
        if not os.path.isfile(path):
            exit(f"'{path}' doesn't exists")

    # if passed all tests return True
    return True


def read_file(path):
    # read the file and return the sequence
    with open(path, 'r') as file:
        seq = file.read()
    # returning the sequence without '\n' at the end
    return seq[:-1]


def analyze_seq(seq):
    # count the number of A, T, C, G and unknown and returns a dictionary

    # making the analyze case insensative
    seq = seq.upper()

    result = {}
    nuc_num = 0

    # counting each nuclaotide and saving its count to result
    for nuc in ["A", "T", "C", "G"]:
        count = seq.count(nuc)
        result[nuc] = count
        nuc_num += count

    # adding to result any unknown charactors
    result["Unknown"] = len(seq) - nuc_num
    result["Total"] = len(seq)

    # return a dict of the countings
    return result


def print_result(result):
    for nuc, count in result.items():
        if nuc == "Total":
            continue
        print(f"{nuc}:{' ' * (8 - len(nuc))} {count} {round(count / result["Total"] * 100, 1)}%")
    
    print(f"Total:{' ' * 3} {result["Total"]}\n\n")


if __name__ == '__main__':
    main()

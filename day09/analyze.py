import re, argparse, os.path
from Bio import SeqIO
from Bio.Seq import Seq



def main():

    # Setting up argparse
    parser = argparse.ArgumentParser(description="Analyze sequence files.")
    parser.add_argument("FILE", help="Path to the input file (FASTA or GenBank).")
    parser.add_argument("--duplicate", action="store_true", help="Look for the longest subsequence that appears twice.")
    parser.add_argument("--reverse_complement", action="store_true", help="Generate the reverse complement of all sequences.")
    
    args = parser.parse_args()

    # check if file exist
    if not os.path.isfile(args.FILE):
        exit(f"'{args.FILE}' doesn't exists")

    # extract the sequences from the file
    sequences = extract_sequences(args.FILE)
    if not sequences:  # validate the extraction went well
        print("No sequences were found in the file.")
        return

    if args.duplicate:      
        for seq_id, seq in sequences.items():
            print(f'The longest subsequence in {seq_id} is: {longest_subseq(seq)}')

    if args.reverse_complement:
        for seq_id, seq in sequences.items():
            reverse_complement = str(Seq(seq).reverse_complement())
            print(f"Reverse complement of {seq_id}: {reverse_complement}")



def extract_sequences(file_path):
    """
    Extracts sequences from a file in FASTA or GenBank format.
    """
    sequences = {}

    try:
        # Use SeqIO to automatically detect format and parse the file
        with open(file_path, "r") as file:
            for record in SeqIO.parse(file, "fasta"):
                sequences[record.id] = str(record.seq)

        if not sequences:  # If no sequences found, try GenBank
            with open(file_path, "r") as file:
                for record in SeqIO.parse(file, "genbank"):
                    sequences[record.id] = str(record.seq)

    except Exception as e:
        print(f"Error: {e}")

    return sequences


def longest_subseq(seq):
    """
    Returns the longest subsequence that appears twice.
    """
    length = 1
    result = ''
    while True:
        regex = r'([GATC]{' + str(length) + r'}).*\1'
        m = re.search(regex, seq)
        if m:
            result = m.group(1)
            length += 1
        else:
            break
    
    return result


if __name__ == "__main__":
    main()

import argparse                              # for extracting CMA
from datetime import datetime                # for logging the date and time of search
from Bio import Entrez                       # to search NCBI
import csv                                   # to write to a csv file the search log
import os                                    # to check if the file already exists or need to initialize it
from urllib.error import HTTPError           # dealing with formating errors

# email for NCBI API
Entrez.email = "rotem.berda@weizmann.ac.il"


def main():

    # setting the cla
    parser = argparse.ArgumentParser()
    parser.add_argument('--database', help= 'Database to search in', required= False, type = str, default= 'nucleotide')
    parser.add_argument('--term', help = 'TERM to search for', required= True, type = str)
    parser.add_argument('--number', help = 'Max results', required= False, type= int, default= 10)

    args = parser.parse_args()

    # search in NCBI 
    search_results = search_ncbi(args.database, args.term, args.number)
    ids = search_results.get("IdList", [])
    total_items = int(search_results.get("Count", 0))

    # fetch and save each file
    file_names = fetch_file_names(ids, args.database, args.term)

    # print file names
    print("Saved files:", ", ".join(file_names))

    # log the search results to a csv
    log_result_to_csv(args.database, args.term, args.number, total_items)

def search_ncbi(database, term, number):
    # search the term in the database and return the results

    handle = Entrez.esearch(db= database, term= term, retmax= number)
    record = Entrez.read(handle)
    return record

def fetch_file_names(ids, database, term):
    # fetching the file names and returning in a list

    # if file path doesn't exist, create new one
    folder_name = "result_files"
    os.makedirs(folder_name, exist_ok=True)

    file_names = []
    for id in ids:
        try:
            handle = Entrez.efetch(db=database, id=id, rettype="fasta", retmode="text")
            file_name = f"{database}_{term}_{id}.fasta"
            file_path = os.path.join(folder_name, file_name)

            # writing the results to a file named 'database_term_id.fasta'
            with open(file_path, "w") as file:
                file.write(handle.read())
            # add the file name to the list
            file_names.append(file_path)

        except HTTPError as e:
            # handles with errors
            print(f"Error fetching ID {id}: {e}")

    return file_names

def log_result_to_csv(database, term, number, total_items):
    # logs the searching to a csv file

    log_file = "ncbi_search_log.csv"

    # getting the date and time of the search
    date_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_data = [date_time, term, database, number, total_items]
    
    if not os.path.isfile(log_file):
        file_exist = False
    else:
        file_exist = True

    with open(log_file, mode="a", newline="") as csvfile:
        writer = csv.writer(csvfile)
        # if file doesn't exist, creating new one with columns names, then adding the search log
        if not file_exist:
            writer.writerow(["date", "term", "database", "max", "total"])
        writer.writerow(log_data)


if __name__ == '__main__':
    main()

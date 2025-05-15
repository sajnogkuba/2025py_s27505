#!/usr/bin/env python3
"""
NCBI GenBank Data Retriever
Podstawowy skrypt do łączenia się z NCBI i pobierania rekordów sekwencji
genetycznych dla danego identyfikatora taksonomicznego.
"""

from Bio import Entrez, SeqIO
import time
import os
import matplotlib.pyplot as plt
import pandas as pd
from io import StringIO

class NCBIRetriever:
    def __init__(self, email, api_key):
        """Initialize with NCBI credentials."""
        self.email = email
        self.api_key = api_key

        # Ustawienia Entrez
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'

    def search_taxid(self, taxid, min_len=None, max_len=None):
        """Search for all records associated with a taxonomic ID and optional length filter."""
        print(f"Searching for records with taxID: {taxid}")
        try:
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"Organism: {organism_name} (TaxID: {taxid})")

            # Budowanie zapytania
            search_term = f"txid{taxid}[Organism]"
            if min_len is not None:
                search_term += f" AND {min_len}:1000000000[Sequence Length]"
            if max_len is not None:
                # jeśli podano oba — zamieniamy zakres
                if min_len is not None:
                    search_term = f"txid{taxid}[Organism] AND {min_len}:{max_len}[Sequence Length]"
                else:
                    search_term += f" AND 0:{max_len}[Sequence Length]"

            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])

            if count == 0:
                print(f"No records found for {organism_name}")
                return None

            print(f"Found {count} records matching length criteria")

            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]
            self.count = count
            return count

        except Exception as e:
            print(f"Error searching TaxID {taxid}: {e}")
            return None

    def fetch_records(self, start=0, max_records=10):
        """Fetch a batch of records using the stored search results."""
        if not hasattr(self, 'webenv') or not hasattr(self,'query_key'):
            print("No search results to fetch. Run search_taxid() first.")
            return []

        try:
            # Limit, aby zapobiec przeciążeniu serwera
            batch_size = min(max_records, 500)

            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )

            # Surowy rekord GenBank
            records_text = handle.read()

            return records_text

        except Exception as e:
            print(f"Error fetching records: {e}")
            return ""

def generate_csv_from_gb_text(gb_text, csv_filename):
    handle = StringIO(gb_text)
    records = SeqIO.parse(handle, "genbank")
    data = [{"accession": r.id, "length": len(r.seq), "description": r.description} for r in records]
    pd.DataFrame(data).to_csv(csv_filename, index=False)
    print(f"Saved {len(data)} records to {csv_filename}")


def generate_plot_from_csv(csv_file, png_file):
    """Generates a line plot of sequence lengths from a CSV file and saves as PNG."""
    df = pd.read_csv(csv_file)
    df = df.sort_values(by="length", ascending=False)

    plt.figure(figsize=(10, 6))
    plt.plot(df["accession"], df["length"], marker='o')
    plt.xticks(rotation=90)
    plt.xlabel("GenBank Accession")
    plt.ylabel("Sequence Length")
    plt.title("Sequence Lengths (sorted)")
    plt.tight_layout()
    plt.savefig(png_file)
    plt.close()
    print(f"Plot saved to {png_file}")


def main():
    # Uzyskaj dane uwierzytelniające
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")

    # Utwórz obiekt retriever
    retriever = NCBIRetriever(email, api_key)

    # Uzyskaj taxid od użytkownika
    taxid = input("Enter taxonomic ID (taxid) of the organism: ")

    # Uzyskaj minimalną i maksymalną długość sekwencji
    min_len = input("Enter minimum sequence length (press Enter to skip): ")
    max_len = input("Enter maximum sequence length (press Enter to skip): ")

    # Konwersja lub ustawienie None
    min_len = int(min_len) if min_len.strip().isdigit() else None
    max_len = int(max_len) if max_len.strip().isdigit() else None

    # Szukaj rekordów z ograniczeniem długości
    count = retriever.search_taxid(taxid, min_len=min_len, max_len=max_len)

    if not count:
        print("No records found. Exiting.")
        return

    # Pobierz kilka pierwszych rekordów jako próbkę
    print("\nFetching sample records...")
    sample_records = retriever.fetch_records(start=0, max_records=50)
    csv_file = f"taxid_{taxid}_sample.csv"
    generate_csv_from_gb_text(sample_records, csv_file)

    png_file = f"taxid_{taxid}_plot.png"
    generate_plot_from_csv(csv_file, png_file)
    print(f"Saved sample records to {csv_file}")
    print("\nNote: This is just a basic retriever. You need to extend its functionality!")
if __name__ == "__main__":
    main()
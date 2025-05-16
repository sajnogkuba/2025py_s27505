from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt
import pandas as pd
from io import StringIO


class Retriever:
    def __init__(self, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'gen'

    def search(self, taxid, min_len=None, max_len=None):
        q = f'txid{taxid}[Organism]'
        if min_len is not None and max_len is not None:
            q = f'txid{taxid}[Organism] AND {min_len}:{max_len}[Sequence Length]'
        elif min_len is not None:
            q += f' AND {min_len}:1000000000[Sequence Length]'
        elif max_len is not None:
            q += f' AND 0:{max_len}[Sequence Length]'
        r = Entrez.read(Entrez.esearch(db='nucleotide', term=q, usehistory='y'))
        if not int(r['Count']): return None
        self.w, self.q, self.c = r['WebEnv'], r['QueryKey'], r['Count']
        return self.c

    def fetch(self, n=5):
        return Entrez.efetch(db='nucleotide', rettype='gb', retmode='text',
                             retstart=0, retmax=n, webenv=self.w, query_key=self.q).read()


def save_csv(text, out):
    d = [{'accession': r.id, 'length': len(r.seq), 'description': r.description}
         for r in SeqIO.parse(StringIO(text), "genbank")]
    pd.DataFrame(d).to_csv(out, index=0)


def save_plot(csv, out):
    d = pd.read_csv(csv).sort_values('length', ascending=False)
    plt.plot(d.accession, d.length, marker='o')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(out)
    plt.close()


def main():
    email = input('email: ')
    key = input('api key: ')
    taxid = input('taxid: ')
    min_len = input('min len: ')
    max_len = input('max len: ')
    min_len = int(min_len) if min_len.isdigit() else None
    max_len = int(max_len) if max_len.isdigit() else None
    r = Retriever(email, key)
    if not r.search(taxid, min_len, max_len): return
    txt = r.fetch()
    csv = f"taxid_{taxid}_sample.csv"
    save_csv(txt, csv)
    save_plot(csv, csv.replace('.csv', '.png'))


main()

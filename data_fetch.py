from Bio import Entrez
from data.pubmed_keywords import KEYWORDS
from tqdm import tqdm
import os
import csv

Entrez.email = "jakub.lechowski@student.uj.edu.pl"


def fetch_abstracts_and_title(keywords, amount, output_file_name):
    print("[DATA FETCH] Fetching...")
    abstracts_and_title = []
    for index, keyword in enumerate(keywords):
        print(f"{index}) {keyword}")
        search_handle = Entrez.esearch(
            db="pubmed", term=keyword, retmax=amount
        )
        search_results = Entrez.read(search_handle)
        ids = set(search_results["IdList"])
        for article_id in tqdm(ids):
            try:
                fetch_handle = Entrez.efetch(db="pubmed", id=article_id, retmode="xml")
                record = Entrez.read(fetch_handle)
                title = record["PubmedArticle"][0]["MedlineCitation"]["Article"][
                    "ArticleTitle"
                ]

                abstract_xml = record["PubmedArticle"][0]["MedlineCitation"][
                    "Article"
                ]["Abstract"]["AbstractText"]
                abstract = "".join([str(x) for x in abstract_xml])
            except:
                continue
            abstracts_and_title.append((abstract, title))
    print(
        f"{len(abstracts_and_title)} Abstract and titles fetched!"
    )

    with open(
        os.path.join("data", f"{output_file_name}.csv"), "w", newline=""
    ) as f:
        writer = csv.writer(f)
        writer.writerows(abstracts_and_title)


fetch_abstracts_and_title(KEYWORDS, 1000, "abstracts_and_titles_10000x100")

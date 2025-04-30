from bs4 import BeautifulSoup
import requests
import time

import numpy as np

def get_pubmed_count(gene):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    time.sleep(1)
    params = {
        "db": "pubmed",
        "term": f'"{gene}"[TIAB] AND cancer[MH]',
        "retmode": "json"  # request JSON response
    }
    response = requests.get(base_url, params=params)
    try:
        response.raise_for_status()
    except Exception as e:
        print('Error Occured')
        return 0
    
    data = response.json()
    # 'count' is the total number of PubMed articles matching the query
    
    # count = int(data["esearchresult"]["count"])
    count = int(data["esearchresult"].get("count", 0))
    
    print(count)
    return count

def get_pubmed_counts(genes):
    
    return {gene: get_pubmed_count(gene) for gene in genes}





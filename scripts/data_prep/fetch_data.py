import os
import requests
import sys
from pathlib import Path

def download_file(url, dest_path):
    if dest_path.exists():
        print(f"File already exists: {dest_path}")
        return

    print(f"Downloading {url} to {dest_path}...")
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(dest_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print("Download complete.")
    except Exception as e:
        print(f"Failed to download {url}: {e}")
        sys.exit(1)

def fetch_datasets():
    # Resolve Project Root (2 levels up from this script: scripts/data_prep/ -> scripts/ -> root)
    script_path = Path(__file__).resolve()
    project_root = script_path.parent.parent.parent
    
    raw_data_dir = project_root / "data" / "raw"
    
    # Ensure directories exist
    (raw_data_dir / "GSE215865").mkdir(parents=True, exist_ok=True)
    (raw_data_dir / "GSE157859").mkdir(parents=True, exist_ok=True)

    print(f"Project Root: {project_root}")
    print(f"Data Directory: {raw_data_dir}")

    # GSE215865 (COVID-19 Blood)
    url1 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE215nnn/GSE215865/suppl/GSE215865_rnaseq_logCPM_matrix.csv.gz"
    dest1 = raw_data_dir / "GSE215865" / "GSE215865_rnaseq_logCPM_matrix.csv.gz"
    download_file(url1, dest1)

    # GSE157859 (Infection Response)
    # Correct filename verified via search: GSE157859_TPM_matrix_of_lncRNA-mRNA_all_Samples.txt.gz
    url2 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157859/suppl/GSE157859_TPM_matrix_of_lncRNA-mRNA_all_Samples.txt.gz"
    dest2 = raw_data_dir / "GSE157859" / "GSE157859_TPM_matrix_of_lncRNA-mRNA_all_Samples.txt.gz"
    download_file(url2, dest2)

if __name__ == "__main__":
    fetch_datasets()

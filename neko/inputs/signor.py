import requests
import pandas as pd
from io import BytesIO

def download_signor_database(save_path: str = None) -> pd.DataFrame:
    """
    Download the SIGNOR human dataset from the official API.
    If save_path is provided, saves the file as TSV. Otherwise, returns a DataFrame.
    """
    url = "https://signor.uniroma2.it/API/getHumanData.php"
    try:
        r = requests.get(url)
        r.raise_for_status()
        if save_path:
            with open(save_path, 'wb') as f:
                f.write(r.content)
            return None
        else:
            df = pd.read_csv(BytesIO(r.content), sep='\t')
            return df
    except requests.RequestException as e:
        raise RuntimeError(f"Error downloading SIGNOR database: {str(e)}")

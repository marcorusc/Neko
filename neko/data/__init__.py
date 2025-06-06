import io
import urllib.parse

import pandas as pd
import requests

BASEURL = "https://commons.omnipathdb.org/"


def _open(url: str, ftype: str | None = None, df: bool | dict = False):
    """Open a remote file and optionally return it as a DataFrame."""

    if not ftype:
        ftype = urllib.parse.urlparse(url).path.split(".")[-1].lower()

    if df is not False:
        df_kwargs = df if isinstance(df, dict) else {}
        if ftype == "tsv":
            return pd.read_table(url, sep="\t", **df_kwargs)
        elif ftype in {"csv", "txt"}:
            return pd.read_csv(url, **df_kwargs)
        elif ftype in {"xls", "xlsx"}:
            return pd.read_excel(url, **df_kwargs)
        else:
            raise NotImplementedError(f"Unsupported file type {ftype}")

    resp = requests.get(url)
    resp.raise_for_status()
    return io.BytesIO(resp.content)


def _baseurl() -> str:
    return BASEURL


def _retrieve(path: str, ftype: str = "tsv") -> pd.DataFrame:

    url = urllib.parse.urljoin(BASEURL, path)

    return _open(url, ftype=ftype, df=True)


def phosphosite_kinase_substrate():

    return _retrieve('phosphosite/kinase-substrate.tsv')


def phosphosite_regulatory_sites():

    return _retrieve('phosphosite/regulatory-sites.tsv')


def huri(dataset: str = 'HI-union', translated = True) -> pd.DataFrame:

    url = (
        'https://github.com/sysbio-curie/Medulloblastoma_project/'
        'raw/main/Huri_analysis/data/%s%s.csv'
    )

    url = url % (dataset, '_translated' if translated else '')

    df = _open(url, df=True)

    return df

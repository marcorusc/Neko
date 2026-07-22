import re
import logging
import os
import uuid

import omnipath as op
import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

import pandas as pd

from neko._cache import cache_dir


_HPA_TISSUE_RESOURCE = "HPA_tissue"
_REQUIRED_ANNOTATION_COLUMNS = frozenset({
    "genesymbol",
    "label",
    "record_id",
    "value",
})
_NOT_DETECTED_LEVELS = frozenset({"not detected"})
_HPA_CANCER_DATA_URL = (
    "https://www.proteinatlas.org/download/tsv/cancer_data.tsv.zip"
)
_HPA_CANCER_CACHE_SUBDIR = "annotations"
_HPA_CANCER_CACHE_FILENAME = "hpa_cancer_data.tsv.zip"
_MIN_HPA_CANCER_ROWS = 100_000
_HPA_CANCER_COLUMNS = frozenset({
    "Gene name",
    "Cancer",
    "High",
    "Medium",
    "Low",
    "Not detected",
})
_HPA_CANCER_TISSUES = frozenset({
    "breast cancer",
    "carcinoid",
    "cervical cancer",
    "colorectal cancer",
    "endometrial cancer",
    "glioma",
    "head and neck cancer",
    "liver cancer",
    "lung cancer",
    "lymphoma",
    "melanoma",
    "ovarian cancer",
    "pancreatic cancer",
    "prostate cancer",
    "renal cancer",
    "skin cancer",
    "stomach cancer",
    "testis cancer",
    "thyroid cancer",
    "urothelial cancer",
})

logger = logging.getLogger(__name__)


class AnnotationServiceError(RuntimeError):
    """Raised when tissue annotations cannot be retrieved or interpreted."""


def _normalize_annotation_value(value):
    """Normalize an OmniPath annotation value for comparison."""
    if pd.isna(value):
        return ""
    return " ".join(str(value).casefold().split())


def _hpa_cancer_cache_path():
    """Return the managed cache path for the HPA cancer table."""
    return (
        cache_dir(_HPA_CANCER_CACHE_SUBDIR)
        / _HPA_CANCER_CACHE_FILENAME
    )


def _new_hpa_session():
    """Create a bounded-retry session for Human Protein Atlas downloads."""
    session = requests.Session()
    retries = Retry(
        total=3,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset({'GET'}),
    )
    session.mount('https://', HTTPAdapter(max_retries=retries))
    return session


def _read_hpa_cancer_table(path):
    """Read and validate the cached HPA cancer expression table."""
    table = pd.read_csv(path, sep='\t', compression='zip', low_memory=False)
    missing_columns = _HPA_CANCER_COLUMNS.difference(table.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"HPA cancer data is missing columns: {missing}.")
    if len(table) < _MIN_HPA_CANCER_ROWS:
        raise ValueError(
            "HPA cancer data is unexpectedly small "
            f"({len(table)} rows)."
        )
    return table


def _download_hpa_cancer_table(timeout=120):
    """Download, validate, and atomically cache the HPA cancer table."""
    path = _hpa_cancer_cache_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    token = f"{os.getpid()}.{uuid.uuid4().hex}"
    temporary_path = path.with_name(f".{path.name}.{token}.tmp")

    try:
        response = _new_hpa_session().get(
            _HPA_CANCER_DATA_URL,
            timeout=(10, timeout),
        )
        response.raise_for_status()
        temporary_path.write_bytes(response.content)
        table = _read_hpa_cancer_table(temporary_path)
        temporary_path.replace(path)
        return table
    finally:
        temporary_path.unlink(missing_ok=True)


def _load_hpa_cancer_table():
    """Load a valid cached HPA table, downloading it only when necessary."""
    path = _hpa_cancer_cache_path()
    if path.is_file():
        try:
            return _read_hpa_cancer_table(path)
        except (OSError, ValueError, pd.errors.ParserError) as exc:
            logger.warning("Ignoring invalid cached HPA cancer data: %s", exc)
    return _download_hpa_cancer_table()


def _hpa_cancer_expressed_genes(gene_symbols, tissue):
    """Return genes detected by HPA IHC in at least one cancer sample."""
    target_tissue = _normalize_annotation_value(tissue)
    table = _load_hpa_cancer_table()
    cancers = table['Cancer'].map(_normalize_annotation_value)
    symbols = table['Gene name'].map(_normalize_annotation_value)
    requested = {
        _normalize_annotation_value(symbol)
        for symbol in gene_symbols
    }
    matches = table[cancers.eq(target_tissue) & symbols.isin(requested)].copy()
    if matches.empty:
        return set()

    detected_counts = matches[['High', 'Medium', 'Low']].apply(
        pd.to_numeric,
        errors='coerce',
    ).fillna(0).sum(axis=1)
    return set(symbols.loc[matches.index[detected_counts.gt(0)]])


def _expressed_genes(annotations_df, tissue):
    """Return normalized symbols expressed in ``tissue`` in HPA records."""
    if not isinstance(annotations_df, pd.DataFrame):
        raise AnnotationServiceError(
            "OmniPath returned an invalid tissue annotation response: "
            f"expected a pandas DataFrame, got {type(annotations_df).__name__}."
        )

    missing_columns = _REQUIRED_ANNOTATION_COLUMNS.difference(
        annotations_df.columns
    )
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise AnnotationServiceError(
            "OmniPath returned an invalid tissue annotation response; "
            f"missing columns: {missing}."
        )

    if annotations_df.empty:
        return set()

    records = annotations_df.copy()
    records["_label"] = records["label"].map(_normalize_annotation_value)
    records = records[records["_label"].isin(("tissue", "level"))]
    if records.empty:
        return set()

    records = (
        records.pivot_table(
            index=["genesymbol", "record_id"],
            columns="_label",
            values="value",
            aggfunc="first",
            observed=True,
        )
        .reset_index()
    )
    if "tissue" not in records or "level" not in records:
        return set()

    target_tissue = _normalize_annotation_value(tissue)
    tissues = records["tissue"].map(_normalize_annotation_value)
    levels = records["level"].map(_normalize_annotation_value)
    detected = levels.ne("") & ~levels.isin(_NOT_DETECTED_LEVELS)
    matching_records = records[tissues.eq(target_tissue) & detected]

    return {
        _normalize_annotation_value(symbol)
        for symbol in matching_records["genesymbol"]
        if not pd.isna(symbol)
    }


def fetch_nodes_from_url(url):

    """
    fetch the nodes in a list from the given geneontology url
    """
    print("Start requesting genes from Gene Ontology")
    print("Fetching from: ", url)
    response = requests.get(url)
    print("Done")
    if not response:
        print("Error fetching the genes, no entry found, please check the id accession code")
        return
    genes = response.text.split("\n")
    genes = genes[:-1]
    genes_unique = list(set(genes))
    return genes_unique


class Ontology:
    """
    class that stores some functionalities to connect phenotypes to nodes
    and to associate information for each node at tissue level

    """
    def __init__(self):
        self.gene_ontology_url = "https://golr-aux.geneontology.io/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=bioentity_label&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=%7C&fq=document_category:%22annotation%22&fq=isa_partof_closure:%22GO:0007049%22&fq=taxon_subset_closure_label:%22Homo%20sapiens%22&fq=type:%22protein%22&fq=annotation_class_label:%22G1/S%20transition%20of%20mitotic%20cell%20cycle%22&facet.field=aspect&facet.field=taxon_subset_closure_label&facet.field=type&facet.field=evidence_subset_closure_label&facet.field=regulates_closure_label&facet.field=isa_partof_closure_label&facet.field=annotation_class_label&facet.field=qualifier&facet.field=annotation_extension_class_closure_label&facet.field=assigned_by&facet.field=panther_family_label&q=*:*"
        self.accession_to_phenotype_dict = {"GO:0010718": "positive regulation of epithelial to mesenchymal transition"}
        return

    def modify_url_ontology(self, new_go_code, new_description):
        """
        Modifies a given URL by replacing a GO code and a descriptive string
        in very specific locations identified by prefixes.

        Parameters:
        - url (str): The original URL to be modified.
        - new_go_code (str): The new GO code to insert into the URL.
        - new_description (str): The new descriptive string to insert into the URL.

        Returns:
        - str: The modified URL.
        """
        # Ensure we work with plain strings
        if not isinstance(new_go_code, str):
            new_go_code = str(new_go_code)
        # Define the prefixes that identify where the replacements should occur,
        go_code_prefix = "isa_partof_closure:%22"
        description_prefix = "annotation_class_label:%22"

        # Patterns to find the exact locations for replacements
        go_code_pattern = re.escape(go_code_prefix) + r"GO:\d{7}"
        description_pattern = re.escape(description_prefix) + r".+?%22"

        # Perform the replacements
        url = re.sub(
            go_code_pattern,
            go_code_prefix + new_go_code,
            self.gene_ontology_url,
            count=1,
        )
        # URL-encode the new description
        new_description_encoded = re.sub(r" ", "%20", new_description)
        url = re.sub(
            description_pattern,
            description_prefix + new_description_encoded + "%22",
            url,
            count=1,
        )

        return url

    def get_markers(self,
                    phenotype: str = None,
                    id_accession: str = None):
        if phenotype and id_accession:
            self.accession_to_phenotype_dict[id_accession] = phenotype
            url = self.modify_url_ontology(id_accession, phenotype)
            genes = fetch_nodes_from_url(url)
            return genes
        elif phenotype and not id_accession:
            matches = [
                acc for acc, description in self.accession_to_phenotype_dict.items()
                if description == phenotype
            ]
            if not matches:
                print("Invalid GO id accession or phenotype description:")
                print("GO if accession used = ", id_accession)
                print("Phenotype description used = ", phenotype)
                return
            id_accession = matches[0]
            url = self.modify_url_ontology(id_accession, phenotype)
            genes = fetch_nodes_from_url(url)
            return genes
        elif id_accession and not phenotype:
            phenotype = self.accession_to_phenotype_dict[id_accession]
            url = self.modify_url_ontology(id_accession, phenotype)
            genes = fetch_nodes_from_url(url)
            return genes
        else:
            print("Invalid GO id accession or phenotype description:")
            print("GO if accession used = ", id_accession)
            print("Phenotype description used = ", phenotype)
            return

    def check_tissue_annotations(self, genes_df, tissue):
        """
        Check whether genes have detected HPA expression in a tissue.

        Args:
        genes_df (DataFrame): DataFrame containing gene symbols.
        tissue (str): Tissue to match exactly after case/whitespace normalization.

        Returns:
        DataFrame: Gene symbols and their detected-expression status.

        Raises:
        AnnotationServiceError: If OmniPath is unavailable or returns an
            unexpected annotation schema.
        """

        if not isinstance(genes_df, pd.DataFrame):
            raise TypeError(
                "genes_df must be a pandas DataFrame containing a "
                "'Genesymbol' column."
            )
        if 'Genesymbol' not in genes_df:
            raise ValueError("genes_df must contain a 'Genesymbol' column.")
        if not isinstance(tissue, str) or not tissue.strip():
            raise ValueError("tissue must be a non-empty string.")

        raw_symbols = genes_df['Genesymbol']
        if raw_symbols.isna().any():
            raise ValueError("genes_df contains a missing gene symbol.")

        gene_symbols = raw_symbols.astype(str).tolist()
        if any(not symbol.strip() for symbol in gene_symbols):
            raise ValueError("genes_df contains an empty gene symbol.")
        if not gene_symbols:
            return pd.DataFrame({
                'Genesymbol': pd.Series(dtype='object'),
                'in_tissue': pd.Series(dtype='bool'),
            })

        unique_symbols = list(dict.fromkeys(gene_symbols))
        normalized_tissue = _normalize_annotation_value(tissue)

        if normalized_tissue in _HPA_CANCER_TISSUES:
            try:
                expressed = _hpa_cancer_expressed_genes(
                    unique_symbols,
                    tissue,
                )
            except Exception as exc:
                raise AnnotationServiceError(
                    "Unable to retrieve cancer expression data from the "
                    "Human Protein Atlas. No genes were classified as absent."
                ) from exc
        else:
            # A single resource-restricted request is faster, cacheable as one
            # unit, and less likely to fail than one request per gene.
            try:
                annotations_df = op.requests.Annotations.get(
                    proteins=unique_symbols,
                    resources=_HPA_TISSUE_RESOURCE,
                )
            except Exception as exc:
                raise AnnotationServiceError(
                    "Unable to retrieve HPA tissue annotations from OmniPath. "
                    "The service may be temporarily unavailable; no genes "
                    "were classified as absent."
                ) from exc
            expressed = _expressed_genes(annotations_df, tissue)

        return pd.DataFrame({
            'Genesymbol': gene_symbols,
            'in_tissue': [
                _normalize_annotation_value(symbol) in expressed
                for symbol in gene_symbols
            ],
        })

import pandas as pd
import pytest

from neko._annotations import gene_ontology
from neko._annotations.gene_ontology import (
    AnnotationServiceError,
    Ontology,
)


def _annotations():
    return pd.DataFrame([
        {
            'genesymbol': 'SRC',
            'record_id': 1,
            'label': 'tissue',
            'value': ' Colon ',
        },
        {
            'genesymbol': 'SRC',
            'record_id': 1,
            'label': 'level',
            'value': 'Medium',
        },
        {
            'genesymbol': 'TP53',
            'record_id': 2,
            'label': 'tissue',
            'value': 'colon',
        },
        {
            'genesymbol': 'TP53',
            'record_id': 2,
            'label': 'level',
            'value': 'Not detected',
        },
        {
            'genesymbol': 'TP53',
            'record_id': 3,
            'label': 'tissue',
            'value': 'liver cancer',
        },
        {
            'genesymbol': 'TP53',
            'record_id': 3,
            'label': 'level',
            'value': 'High',
        },
    ])


def test_non_cancer_annotations_use_one_omnipath_request(monkeypatch):
    calls = []

    def get(**kwargs):
        calls.append(kwargs)
        return _annotations()

    monkeypatch.setattr(
        gene_ontology.op.requests.Annotations,
        'get',
        get,
    )
    genes = pd.DataFrame({
        'Genesymbol': ['SRC', 'TP53', 'UNKNOWN', 'SRC'],
    })

    result = Ontology().check_tissue_annotations(
        genes,
        '  COLON ',
    )

    assert calls == [{
        'proteins': ['SRC', 'TP53', 'UNKNOWN'],
        'resources': 'HPA_tissue',
    }]
    assert result.to_dict('records') == [
        {'Genesymbol': 'SRC', 'in_tissue': True},
        {'Genesymbol': 'TP53', 'in_tissue': False},
        {'Genesymbol': 'UNKNOWN', 'in_tissue': False},
        {'Genesymbol': 'SRC', 'in_tissue': True},
    ]


def test_tissue_match_is_not_a_substring(monkeypatch):
    monkeypatch.setattr(
        gene_ontology.op.requests.Annotations,
        'get',
        lambda **_: _annotations(),
    )

    result = Ontology().check_tissue_annotations(
        pd.DataFrame({'Genesymbol': ['SRC']}),
        'col',
    )

    assert result['in_tissue'].tolist() == [False]


def test_cancer_annotations_use_direct_hpa_data(monkeypatch):
    table = pd.DataFrame([
        {
            'Gene name': 'SRC',
            'Cancer': 'colorectal cancer',
            'High': 0,
            'Medium': 2,
            'Low': 0,
            'Not detected': 1,
        },
        {
            'Gene name': 'TP53',
            'Cancer': 'colorectal cancer',
            'High': 0,
            'Medium': 0,
            'Low': 0,
            'Not detected': 3,
        },
    ])
    monkeypatch.setattr(
        gene_ontology,
        '_load_hpa_cancer_table',
        lambda: table,
    )
    monkeypatch.setattr(
        gene_ontology.op.requests.Annotations,
        'get',
        lambda **_: pytest.fail('Cancer data must bypass OmniPath'),
    )

    result = Ontology().check_tissue_annotations(
        pd.DataFrame({'Genesymbol': ['SRC', 'TP53', 'UNKNOWN']}),
        ' Colorectal  Cancer ',
    )

    assert result.to_dict('records') == [
        {'Genesymbol': 'SRC', 'in_tissue': True},
        {'Genesymbol': 'TP53', 'in_tissue': False},
        {'Genesymbol': 'UNKNOWN', 'in_tissue': False},
    ]


def test_hpa_cancer_loader_reuses_valid_cache(tmp_path, monkeypatch):
    monkeypatch.setenv('NEKO_CACHE_DIR', str(tmp_path / 'cache'))
    monkeypatch.setattr(gene_ontology, '_MIN_HPA_CANCER_ROWS', 1)
    path = gene_ontology._hpa_cancer_cache_path()
    path.parent.mkdir(parents=True)
    pd.DataFrame([{
        'Gene name': 'SRC',
        'Cancer': 'colorectal cancer',
        'High': 0,
        'Medium': 1,
        'Low': 0,
        'Not detected': 0,
    }]).to_csv(
        path,
        sep='\t',
        index=False,
        compression={
            'method': 'zip',
            'archive_name': 'cancer_data.tsv',
        },
    )
    monkeypatch.setattr(
        gene_ontology,
        '_download_hpa_cancer_table',
        lambda: pytest.fail('A valid cache must not be downloaded again'),
    )

    table = gene_ontology._load_hpa_cancer_table()

    assert table['Gene name'].tolist() == ['SRC']


def test_tissue_annotations_wraps_download_failure(monkeypatch):
    def fail(**_):
        raise RuntimeError('No active exception to reraise')

    monkeypatch.setattr(
        gene_ontology.op.requests.Annotations,
        'get',
        fail,
    )

    with pytest.raises(AnnotationServiceError, match='temporarily unavailable') as error:
        Ontology().check_tissue_annotations(
            pd.DataFrame({'Genesymbol': ['SRC']}),
            'colon',
        )

    assert isinstance(error.value.__cause__, RuntimeError)


def test_tissue_annotations_rejects_invalid_response(monkeypatch):
    monkeypatch.setattr(
        gene_ontology.op.requests.Annotations,
        'get',
        lambda **_: pd.DataFrame({'genesymbol': ['SRC']}),
    )

    with pytest.raises(AnnotationServiceError, match='missing columns'):
        Ontology().check_tissue_annotations(
            pd.DataFrame({'Genesymbol': ['SRC']}),
            'colon',
        )


def test_tissue_annotations_accepts_empty_gene_dataframe(monkeypatch):
    def fail(**_):
        pytest.fail('An empty input must not call OmniPath')

    monkeypatch.setattr(
        gene_ontology.op.requests.Annotations,
        'get',
        fail,
    )

    result = Ontology().check_tissue_annotations(
        pd.DataFrame({'Genesymbol': []}),
        'colorectal cancer',
    )

    assert result.empty
    assert result.dtypes.to_dict() == {
        'Genesymbol': object,
        'in_tissue': bool,
    }


@pytest.mark.parametrize(
    ('genes', 'error'),
    [
        (pd.DataFrame(), "must contain a 'Genesymbol' column"),
        (pd.DataFrame({'Genesymbol': [None]}), 'missing gene symbol'),
        (pd.DataFrame({'Genesymbol': ['  ']}), 'empty gene symbol'),
    ],
)
def test_tissue_annotations_validates_gene_symbols(genes, error):
    with pytest.raises(ValueError, match=error):
        Ontology().check_tissue_annotations(genes, 'colorectal cancer')


@pytest.mark.parametrize('tissue', [None, '', '  '])
def test_tissue_annotations_validates_tissue(tissue):
    with pytest.raises(ValueError, match='non-empty string'):
        Ontology().check_tissue_annotations(
            pd.DataFrame({'Genesymbol': ['SRC']}),
            tissue,
        )

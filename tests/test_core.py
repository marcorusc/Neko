import pandas as pd
import pytest

from neko.inputs import identifier_mapping


def identity(x):
    return x


def setup_module(module):
    module.orig_to_uniprot = identifier_mapping.to_uniprot
    module.orig_to_genesymbol = identifier_mapping.to_genesymbol
    identifier_mapping.to_uniprot = identity
    identifier_mapping.to_genesymbol = identity


def teardown_module(module):
    identifier_mapping.to_uniprot = module.orig_to_uniprot
    identifier_mapping.to_genesymbol = module.orig_to_genesymbol


def sample_df():
    return pd.DataFrame({
        'source': ['P1', 'P2', 'P3'],
        'target': ['P2', 'P3', 'P1'],
        'is_directed': [True, True, True],
        'is_stimulation': [True, False, True],
        'is_inhibition': [False, True, False],
        'form_complex': [False, False, False],
    })


def test_universe_basic():
    from neko.inputs._universe import Universe
    df = sample_df()
    u = Universe(df)
    assert len(u) == 3
    assert {'P1', 'P2', 'P3'} <= u.nodes


def test_network_add_edge():
    from neko.core.network import Network
    df = sample_df()
    net = Network(initial_nodes=['P1'], resources=df)
    edge = pd.DataFrame({
        'source': ['P1'],
        'target': ['P2'],
        'type': ['activation'],
        'references': ['ref'],
        'is_stimulation': [True],
        'is_inhibition': [False],
    })
    net.add_edge(edge)
    assert len(net.edges) == 1
    assert {'P1', 'P2'} <= set(net.edges[['source', 'target']].stack())


def test_network_adds_typed_signor_entity_without_translation():
    from neko.core.tools import mapping_node_identifier
    from neko.core.network import Network

    entity = 'PHENOTYPE:Cell_death'
    resources = sample_df()
    resources.loc[len(resources)] = {
        'source': 'P1',
        'target': entity,
        'is_directed': True,
        'is_stimulation': True,
        'is_inhibition': False,
        'form_complex': False,
    }
    net = Network(initial_nodes=['P1'], resources=resources)

    assert mapping_node_identifier(entity) == [None, entity, entity]
    assert net.add_node(entity)
    matches = net.nodes[['Genesymbol', 'Uniprot']].eq([entity, entity])
    assert matches.all(axis=1).any()


def test_visualizer_wraps_typed_signor_entities():
    from neko._visual.visualize_network import wrap_node_name

    assert wrap_node_name('COMPLEX:P23511_Q13952') == 'P23511_Q13952'
    assert wrap_node_name('COMPLEX_NAME:New complex') == (
        'COMPLEX_NAME_New complex'
    )
    assert wrap_node_name('PHENOTYPE:Cell_death') == 'PHENOTYPE_Cell_death'


def test_exports(tmp_path):
    from neko.core.network import Network
    from neko._outputs.exports import Exports
    df = sample_df()
    net = Network(initial_nodes=['P1'], resources=df)
    edge = pd.DataFrame({
        'source': ['P1'],
        'target': ['P2'],
        'type': ['activation'],
        'references': ['ref'],
        'is_stimulation': [True],
        'is_inhibition': [False],
    })
    net.add_edge(edge)
    exp = Exports(net)
    sif_file = tmp_path / 'test.sif'
    bnet_file = tmp_path / 'test.bnet'
    exp.export_sif(str(sif_file))
    exp.export_bnet(str(bnet_file))
    assert sif_file.exists()
    created = list(tmp_path.glob('test*.bnet'))
    assert created


"""Tests for the Universe class improvements."""
import pytest
import pandas as pd


def sample_df():
    """Create a sample DataFrame for testing."""
    return pd.DataFrame({
        'source': ['P1', 'P2', 'P3'],
        'target': ['P2', 'P3', 'P1'],
        'is_directed': [True, True, True],
        'is_stimulation': [True, False, True],
        'is_inhibition': [False, True, False],
        'form_complex': [False, False, False],
    })


def sample_df2():
    """Create another sample DataFrame for testing."""
    return pd.DataFrame({
        'source': ['X', 'Y'],
        'target': ['Y', 'Z'],
        'is_directed': [True, True],
        'is_stimulation': [True, True],
        'is_inhibition': [False, False],
        'form_complex': [False, False],
    })


# Skip tests that require pypath imports if network is unavailable
try:
    from neko.inputs._universe import Universe
    PYPATH_AVAILABLE = True
except (TypeError, ImportError, OSError):
    PYPATH_AVAILABLE = False


@pytest.mark.skipif(not PYPATH_AVAILABLE, reason="pypath import failed")
class TestUniverse:
    """Test suite for Universe class."""

    def test_create_with_dataframe(self):
        """Test creating Universe with a DataFrame."""
        df = sample_df()
        u = Universe(df, name='test1')
        assert len(u) == 3
        assert 'test1' in u.list_resources()

    def test_list_resources(self):
        """Test listing resources."""
        df = sample_df()
        u = Universe(df, name='test1')
        resources = u.list_resources()
        assert isinstance(resources, list)
        assert 'test1' in resources

    def test_add_resources(self):
        """Test adding resources."""
        df1 = sample_df()
        df2 = sample_df2()
        u = Universe(df1, name='test1')
        u.add_resources(df2, name='test2')
        u.build()
        assert 'test1' in u.list_resources()
        assert 'test2' in u.list_resources()
        assert len(u) == 5

    def test_remove_resources(self):
        """Test removing resources."""
        df1 = sample_df()
        df2 = sample_df2()
        u = Universe(df1, name='test1')
        u.add_resources(df2, name='test2')
        u.build()

        result = u.remove_resources('test2')
        assert result is True
        assert 'test2' not in u.list_resources()
        assert len(u) == 3

    def test_remove_nonexistent_resource(self):
        """Test removing a resource that doesn't exist."""
        df = sample_df()
        u = Universe(df, name='test1')
        result = u.remove_resources('nonexistent')
        assert result is False

    def test_get_resource(self):
        """Test getting a specific resource."""
        df = sample_df()
        u = Universe(df, name='test1')
        resource = u.get_resource('test1')
        assert resource is not None
        assert isinstance(resource, pd.DataFrame)

    def test_get_nonexistent_resource(self):
        """Test getting a resource that doesn't exist."""
        df = sample_df()
        u = Universe(df, name='test1')
        resource = u.get_resource('nonexistent')
        assert resource is None

    def test_clear_resources(self):
        """Test clearing all resources."""
        df = sample_df()
        u = Universe(df, name='test1')
        u.clear_resources()
        assert len(u.list_resources()) == 0
        assert len(u) == 0

    def test_create_with_none(self):
        """Test creating Universe with None."""
        u = Universe(None)
        assert len(u) == 0

    def test_column_mapping(self):
        """Test column mapping when adding resources."""
        df_renamed = pd.DataFrame({
            'src': ['P', 'Q'],
            'tgt': ['Q', 'R'],
            'directed': [True, True],
            'stim': [True, False],
            'inhib': [False, True],
            'cplx': [False, False],
        })
        u = Universe(df_renamed, name='mapped', columns={
            'src': 'source',
            'tgt': 'target',
            'directed': 'is_directed',
            'stim': 'is_stimulation',
            'inhib': 'is_inhibition',
            'cplx': 'form_complex',
        })
        assert 'source' in u.interactions.columns
        assert 'target' in u.interactions.columns

    def test_nodes_property(self):
        """Test the nodes property."""
        df = sample_df()
        u = Universe(df, name='test1')
        nodes = u.nodes
        assert isinstance(nodes, set)
        assert {'P1', 'P2', 'P3'} <= nodes

    def test_contains(self):
        """Test the __contains__ method."""
        df = sample_df()
        u = Universe(df, name='test1')
        assert 'P1' in u
        assert 'NonExistent' not in u

    def test_repr(self):
        """Test the __repr__ method."""
        df = sample_df()
        u = Universe(df, name='test1')
        repr_str = repr(u)
        assert 'Universe' in repr_str
        assert 'test1' in repr_str

    def test_add_resources_with_dict(self):
        """Test adding resources as a dictionary."""
        df1 = sample_df()
        df2 = sample_df2()
        resources_dict = {'resource1': df1, 'resource2': df2}
        u = Universe(resources_dict)
        assert 'resource1' in u.list_resources()
        assert 'resource2' in u.list_resources()

    def test_undirected_resources(self):
        """Test adding undirected resources."""
        df = sample_df()
        u = Universe(df, name='undirected', directed=False)
        # Undirected edges should create mutual edges
        # The original 3 edges should become 6 edges (both directions)
        # minus duplicates
        assert len(u) > 3

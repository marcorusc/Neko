from graphviz import Digraph
from IPython.display import display
from yfiles_jupyter_graphs import GraphWidget
from typing import Dict


def wrap_node_name(node_name):
    if ":" in node_name:
        node_name = node_name.replace(":", "_")
    if node_name.startswith("COMPLEX"):
        return node_name[8:]
    else:
        return node_name


class NetworkVisualizer:
    def __init__(self, network, predefined_node=None, color_by="Effect", noi=False):
        net = network.copy()
        self.dataframe_edges = net.convert_edgelist_into_genesymbol().copy()
        self.dataframe_nodes = net.nodes
        self.initial_nodes = net.initial_nodes
        self.color_by = color_by
        self.noi = noi  # nodes of interest
        self.graph = Digraph(format='pdf')
        self.edge_colors = {
            'stimulation': 'green',
            'inhibition': 'red',
            'form complex': 'blue',
            'bimodal': 'purple',  # Adding the mapping for "bimodal" effect
            'undefined': 'gray',  # Adding the mapping for "undefined" effect
            # Add more custom mappings if needed
        }
        self.node_colors = {}  # Dictionary to store custom node colors
        self.predefined_node = wrap_node_name(predefined_node) if predefined_node else None
        # Apply wrap_node_name function to node names in dataframe_nodes
        self.dataframe_nodes['Genesymbol'] = self.dataframe_nodes['Genesymbol'].apply(wrap_node_name)
        self.dataframe_nodes['Uniprot'] = self.dataframe_nodes['Uniprot'].apply(wrap_node_name)

        # Apply wrap_node_name function to node names in "source" and "target" columns of dataframe_edges
        self.dataframe_edges['source'] = self.dataframe_edges['source'].apply(wrap_node_name)
        self.dataframe_edges['target'] = self.dataframe_edges['target'].apply(wrap_node_name)

        self.add_edges_to_graph()
        self.add_nodes_to_graph()

    def set_custom_edge_colors(self, custom_edge_colors):
        # Update the edge_colors dictionary with custom mappings
        self.edge_colors.update(custom_edge_colors)

    def set_node_colors(self, node_colors):
        # Update the node_colors dictionary with custom node colorsdataframe_nodes
        self.node_colors.update(node_colors)

    def add_edges_to_graph(self):
        for _, row in self.dataframe_edges.iterrows():
            effect = row['Effect']
            source = wrap_node_name(row['source'])
            target = wrap_node_name(row['target'])

            # Display only edges connected to the predefined node
            if self.predefined_node and (source != self.predefined_node and target != self.predefined_node):
                continue

            # Determine edge attributes based on effect
            if effect == 'stimulation':
                arrowhead = 'normal'
                color = 'green'
                dir = 'forward'
            elif effect == 'inhibition':
                arrowhead = 'tee'
                color = 'red'
                dir = 'forward'
            elif effect == 'form complex':
                arrowhead = 'dot'
                color = 'blue'
                dir = 'forward'
            elif effect == 'bimodal':
                arrowhead = 'diamond'
                color = 'purple'
                dir = 'forward'
            else:
                arrowhead = 'normal'  # Default arrow
                color = 'black'
                dir = 'none'

            # Add the edge to the graph with specified attributes
            self.graph.edge(source, target, color=color, arrowhead=arrowhead, dir=dir)

    def add_nodes_to_graph(self):
        for _, row in self.dataframe_nodes.iterrows():
            node = row['Genesymbol']
            # add function to set color
            node_color = 'lightgray'
            if node in self.initial_nodes and self.noi:
                node_color = 'lightyellow'
            node_color = self.node_colors.get(node, node_color)

            # Display only the predefined node and its connections
            if self.predefined_node and (node != self.predefined_node):
                continue

            wrapped_node = wrap_node_name(node)
            self.graph.node(wrapped_node, style='filled', fillcolor=node_color)

    def tissue_mapping(self, tissue_df):
        """
        Color the nodes based on their expression in the tissue of interest (based on data from The Human Protein Atlas).

        Args:
            tissue_df (DataFrame): DataFrame containing results indicating whether each gene symbol has tissue annotations containing the selected tissue.
        """
        for _, row in tissue_df.iterrows():
            gene_symbol = row['Genesymbol']
            in_tissue = row['in_tissue']
            node_color = 'lightblue' if in_tissue else 'lightgray'
            self.node_colors[gene_symbol] = node_color

    def render(self, output_file='network', view=False, highlight_nodes=None, highlight_color='lightyellow'):
        """
        Render the graph.

        Args:
            output_file (str): The name of the output file.
            view (bool): If True, display the graph.
            highlight_nodes (list): List of nodes to highlight.
            highlight_color (str): Color to use for highlighting nodes.
        """
        # If highlight_nodes is provided, set the color for each node in the list
        if highlight_nodes is not None:
            for node in highlight_nodes:
                # first check that the node is in the node dataframe
                if wrap_node_name(node) in self.dataframe_nodes['Genesymbol'].values:
                    # then change the color only of the node in the graph
                    self.graph.node(wrap_node_name(node), style='filled', fillcolor=highlight_color)
        if view:
            self.graph.view(filename=output_file)
        else:
            self.graph.render(filename=output_file)
            display(self.graph)  # Display the graph directly in the Jupyter Notebook

    def yfiles_visual(
        self,
        graph_layout,
        directed,
    ):
        # creating empty object for visualization
        w = GraphWidget()

        # filling w with nodes
        objects = []
        for idx, item in self.dataframe_nodes.iterrows():
            obj = {
                "id": self.dataframe_nodes["Uniprot"].loc[idx],
                "properties": {"label": self.dataframe_nodes["Genesymbol"].loc[idx]},
                "color": "#ffffff",
                "styles": {"backgroundColor": "#ffffff"}
            }
            objects.append(obj)
        w.nodes = objects

        # filling w with edges
        objects = []
        for index, row in self.dataframe_edges.iterrows():
            obj = {
                "id": self.dataframe_edges["Effect"].loc[index],
                "start": self.dataframe_edges["source"].loc[index],
                "end": self.dataframe_edges["target"].loc[index],
                "properties": {"references": self.dataframe_edges["References"].loc[index]}}
            objects.append(obj)
        w.edges = objects

        def custom_edge_color_mapping(edge: Dict):
            """let the edge be red if the interaction is an inhibition, else green"""
            return ("#fa1505" if edge['id'] == "inhibition" else "#05e60c")

        w.set_edge_color_mapping(custom_edge_color_mapping)

        def custom_node_color_mapping(node: Dict):
            return {"color": "#ffffff"}

        w.set_node_styles_mapping(custom_node_color_mapping)

        def custom_factor_mapping(node: Dict):
            """choose random factor"""
            return 5

        w.set_node_scale_factor_mapping(custom_factor_mapping)

        def custom_label_styles_mapping(node: Dict):
            """let the label be the negated purple big index"""
            return {
                'text': node["properties"]["label"],
                'backgroundColor': None,
                'fontSize': 40,
                'color': '#030200',
                'shape': 'round-rectangle',
                'textAlignment': 'center'
            }

        w.set_node_label_mapping(custom_label_styles_mapping)

        w.directed = directed
        w.graph_layout = graph_layout

        display(w)

    def vis_comparison(
        self,
        int_comparison,
        node_comparison,
        graph_layout,
        directed,
    ):
        # creating empty object for visualization
        w = GraphWidget()

        objects = []
        for idx, item in node_comparison.iterrows():
            obj = {
                "id": node_comparison["node"].loc[idx],
                "properties": {"label": node_comparison["node"].loc[idx],
                               "comparison": node_comparison["comparison"].loc[idx], },
                "color": "#ffffff",
                #       "styles":{"backgroundColor":"#ffffff"}
            }
            objects.append(obj)
        w.nodes = objects

        # filling w with edges
        objects = []
        for index, row in int_comparison.iterrows():
            obj = {
                "id": int_comparison["comparison"].loc[index],
                "properties": {
                    "comparison": int_comparison["comparison"].loc[index]},
                "start": int_comparison["source"].loc[index],
                "end": int_comparison["target"].loc[index]
            }
            objects.append(obj)
        w.edges = objects

        def custom_node_color_mapping(node: Dict):
            if node['properties']['comparison'] == "Unique to Network 1":
                return {"color": "#f5f536"}
            elif node['properties']['comparison'] == "Unique to Network 2":
                return {"color": "#36f55f"}
            elif node['properties']['comparison'] == "Common":
                return {"color": "#3643f5"}

        w.set_node_styles_mapping(custom_node_color_mapping)

        def custom_factor_mapping(node: Dict):
            """choose random factor"""
            return 5

        w.set_node_scale_factor_mapping(custom_factor_mapping)

        def custom_label_styles_mapping(node: Dict):
            """let the label be the negated purple big index"""
            return {
                'text': node["id"],
                'backgroundColor': None,
                'fontSize': 20,
                'color': '#030200',
                'position': 'center',
                'maximumWidth': 130,
                'wrapping': 'word',
                'textAlignment': 'center'
            }

        w.set_node_label_mapping(custom_label_styles_mapping)

        def custom_edge_color_mapping(edge: Dict):
            if edge['id'] == "Unique to Network 1":
                return "#e3941e"
            elif edge['id'] == "Unique to Network 2":
                return "#36f55f"
            elif edge['id'] == "Common":
                return "#3643f5"
            elif edge['id'] == "Conflicting":
                return "#ffcc00"

        w.set_edge_color_mapping(custom_edge_color_mapping)

        w.directed = directed
        w.graph_layout = graph_layout

        display(w)

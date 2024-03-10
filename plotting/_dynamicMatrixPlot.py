from scanpy.plotting._matrixplot import MatrixPlot
import plotly.graph_objects as go
import scanpy as sc
import pandas as pd
import numpy as np
import warnings
import ipywidgets as widgets
from IPython.display import display

class DynamicMatrixPlot(MatrixPlot):
    def __init__(self, adata, var_names, groupby, **kwds):
        super().__init__(adata, var_names, groupby, **kwds)
        self.output = widgets.Output()

    def compute_stats(self, marker_genes_dict, pbmc):
        # Calcul des moyennes d'expression des gènes marqueurs
        mean_expressions = []
        for group, genes in marker_genes_dict.items():
            expression_data = sc.get.obs_df(pbmc, keys=genes + ['clusters'])
            mean_expression = expression_data.groupby('clusters')[genes].mean()
            mean_expression.columns = [f"{group}_{gene}" for gene in genes]
            mean_expressions.append(mean_expression)

        self.mean_expression_df = pd.concat(mean_expressions, axis=1)

        # Calcul des statistiques supplémentaires
        stats_expressions = []
        for group, genes in marker_genes_dict.items():
            expression_data = sc.get.obs_df(pbmc, keys=genes + ['clusters'])
            stats = expression_data.groupby('clusters')[genes].agg(['mean', 'median', 'std'])
            stats.columns = [f"{group}_{stat}_{gene}" for gene in genes for stat in ['mean', 'median', 'std']]
            stats_expressions.append(stats)

        self.stats_expression_df = pd.concat(stats_expressions, axis=1)

        # Préparation des tooltips pour chaque cellule
        tooltip_texts = self.mean_expression_df.copy()
        for col in self.mean_expression_df.columns:
            group, gene = col.split('_')[0], col.split('_')[-1]
            mean_col = f"{group}_mean_{gene}"
            median_col = f"{group}_median_{gene}"
            std_col = f"{group}_std_{gene}"
            
            for index in self.mean_expression_df.index:
                mean = self.stats_expression_df.at[index, mean_col]
                median = self.stats_expression_df.at[index, median_col]
                std = self.stats_expression_df.at[index, std_col]
                tooltip_texts.at[index, col] = (f"<b>{gene}</b><br>"
                                                f"Mean: {mean:.2f}<br>"
                                                f"Median: {median:.2f}<br>"
                                                f"Std: {std:.2f}")

        self.tooltip_texts = tooltip_texts

    def display_widgets(self, marker_genes_dict):
        self.group_selector = widgets.SelectMultiple(options=list(marker_genes_dict.keys()), description='Groups:', rows=10)
        display_button = widgets.Button(description='Update Graph')
        display_button.on_click(self.update_plot)
        display(widgets.VBox([widgets.HBox([self.group_selector, display_button]), self.output]))

    def update_plot(self, button):
        selected_groups = self.group_selector.value
        self.plot_with_plotly(selected_groups)

    def plot_with_plotly(self, selected_groups):
        with self.output:
            self.output.clear_output(wait=True)  # Efface l'ancienne sortie
            if 'All' in selected_groups or len(selected_groups) == 0:
                filtered_df = self.mean_expression_df
                filtered_tooltips = self.tooltip_texts
            else:
                filtered_cols = [col for col in self.mean_expression_df.columns if col.split('_')[0] in selected_groups]
                filtered_df = self.mean_expression_df[filtered_cols]
                filtered_tooltips = self.tooltip_texts[filtered_cols]

            fig = go.Figure(data=go.Heatmap(
                z=filtered_df.values,
                x=filtered_df.columns,
                y=filtered_df.index,
                text=filtered_tooltips.values,
                hoverongaps=False,
                colorscale='magma'))

            fig.update_traces(hoverinfo='text')
            fig.update_layout(
                title='Average Expression of Marker Genes per Cluster',
                xaxis_title='Marker Gene',
                yaxis_title='Cluster',
                hovermode='closest',
                xaxis={'side': 'bottom'},
                yaxis=dict(tickmode='array', tickvals=np.arange(len(filtered_df.index))),
                template='plotly_white')

            fig.show()

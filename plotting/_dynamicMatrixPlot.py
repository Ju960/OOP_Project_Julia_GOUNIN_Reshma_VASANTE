# importation 
import plotly.graph_objects as go
import scanpy as sc
from scanpy.plotting._matrixplot import MatrixPlot
import pandas as pd
import numpy as np

#Class for user interface
class DynamicMatrixPlot(MatrixPlot):
    def __init__(self, adata, var_names, groupby, **kwargs):
        super().__init__(adata, var_names, groupby, **kwargs)
        self.adata = adata
        self.var_names = var_names

    def plotly_heatmap(self, marker_genes_dict):
        # Calcul of mean expression of marker genes and statistics
        mean_expressions = []
        stats_expressions = []
        for group, genes in marker_genes_dict.items():
            expression_data = sc.get.obs_df(self.adata, keys=genes + ['clusters'])
            mean_expression = expression_data.groupby('clusters')[genes].mean()
            mean_expression.columns = [f"{group}_{gene}" for gene in genes]
            mean_expressions.append(mean_expression)

            stats = expression_data.groupby('clusters')[genes].agg(['mean', 'median', 'std'])
            stats.columns = [f"{group}_{stat}_{gene}" for gene in genes for stat in ['mean', 'median', 'std']]
            stats_expressions.append(stats)

        mean_expression_df = pd.concat(mean_expressions, axis=1)
        stats_expression_df = pd.concat(stats_expressions, axis=1)

        # Preparation of tooltips
        tooltip_texts = mean_expression_df.copy()
        for col in mean_expression_df.columns:
            group, gene = col.split('_')[0], col.split('_')[-1]
            mean_col = f"{group}_mean_{gene}"
            median_col = f"{group}_median_{gene}"
            std_col = f"{group}_std_{gene}"
            
            for index in mean_expression_df.index:
                mean = stats_expression_df.at[index, mean_col]
                median = stats_expression_df.at[index, median_col]
                std = stats_expression_df.at[index, std_col]
                tooltip_texts.at[index, col] = (f"<b>{gene}</b><br>"
                                                f"Mean: {mean:.2f}<br>"
                                                f"Median: {median:.2f}<br>"
                                                f"Std: {std:.2f}")

        # Creation of the Plotly heatmap
        fig = go.Figure(data=go.Heatmap(
                z=mean_expression_df.values,
                x=mean_expression_df.columns,
                y=mean_expression_df.index,
                text=tooltip_texts.values,
                hoverongaps=False,
                colorscale='magma'))

        fig.update_traces(hoverinfo='text')
        fig.update_layout(
            title='Average Expression of Marker Genes per Cluster',
            xaxis_title='Marker Gene',
            yaxis_title='Cluster',
            hovermode='closest',
            xaxis={'side': 'bottom'},
            yaxis=dict(tickmode='array', tickvals=np.arange(len(mean_expression_df.index))),
            template='plotly_white')

        return fig


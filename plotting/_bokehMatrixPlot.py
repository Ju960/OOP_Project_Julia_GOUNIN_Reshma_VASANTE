import numpy as np
import pandas as pd
import scanpy as sc
from bokeh.models import HoverTool, ColorBar, ColumnDataSource
from bokeh.plotting import figure, output_notebook, show
from bokeh.transform import linear_cmap
from scanpy.plotting._matrixplot import MatrixPlot

# Enabling Bokeh display in the notebook
output_notebook()

class BokehMatrixPlot(MatrixPlot):
    def __init__(self, adata, var_names, groupby, use_raw=None, **kwargs):
        # Ensure that var_names is a list
        if not isinstance(var_names, list):
            raise ValueError("var_names doit être une liste")
        super().__init__(adata, var_names, groupby, use_raw=use_raw, **kwargs)
    
    def to_bokeh(self):
        # Preparation of data
        data_frames = []
        for gene in self.var_names:
            expression_data = sc.get.obs_df(self.adata, keys=[gene] + ['clusters'])
            temp_df = expression_data[['clusters', gene]]
            temp_df['gene'] = gene
            temp_df.rename(columns={gene: 'expression'}, inplace=True)
            data_frames.append(temp_df)

        df_bokeh = pd.concat(data_frames)
        grouped_df = df_bokeh.groupby(['clusters', 'gene'], as_index=False).agg(
            expression_mean=('expression', 'mean'),
            expression_min=('expression', 'min'),
            expression_max=('expression', 'max'),
        )
        
        unique_clusters = sorted(grouped_df['clusters'].unique().tolist())
        unique_genes = sorted(grouped_df['gene'].unique().tolist())

        source = ColumnDataSource(grouped_df)
        mapper = linear_cmap(field_name='expression_mean', palette='Blues256', low=grouped_df.expression_mean.min(), high=grouped_df.expression_mean.max())

        p = figure(x_range=unique_clusters, y_range=unique_genes, title="Expression of Genes per Cluster",
                   x_axis_label='Cluster', y_axis_label='Genes', tools="")

        tooltips = """
        <div>
            <span style="font-size: 12px; font-weight: bold;">Cluster:</span>
            <span style="font-size: 12px;">@clusters</span>
            <span style="margin-left: 15px; font-size: 12px; font-weight: bold;">Gene:</span>
            <span style="font-size: 12px;">@gene</span>
            <span style="margin-left: 15px; font-size: 12px; font-weight: bold;">Mean Expression:</span>
            <span style="font-size: 12px;">@expression_mean{0.2f}</span>
            <span style="margin-left: 15px; font-size: 12px; font-weight: bold;">Min Expression:</span>
            <span style="font-size: 12px;">@expression_min{0.2f}</span>
            <span style="margin-left: 15px; font-size: 12px; font-weight: bold;">Max Expression:</span>
            <span style="font-size: 12px;">@expression_max{0.2f}</span>
        </div>
        """
        
        p.rect(x='clusters', y='gene', width=1, height=1, source=source, line_color=None, fill_color=mapper)
        
        color_bar = ColorBar(color_mapper=mapper['transform'], width=8, location=(0,0))
        p.add_layout(color_bar, 'right')
        
        p.add_tools(HoverTool(tooltips=tooltips))
        
        show(p)
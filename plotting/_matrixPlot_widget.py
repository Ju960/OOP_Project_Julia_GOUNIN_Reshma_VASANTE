import scanpy as sc
import ipywidgets as widgets
from IPython.display import display
from scanpy.plotting._matrixplot import MatrixPlot


# Class for User interface
class InteractiveMatrixPlot(MatrixPlot):
    def __init__(self, adata, var_names, **kwargs):
        super().__init__(adata, var_names, groupby='cell type', **kwargs)
        self.adata = adata
        self.var_names = var_names

        # Creation of widgets
        self.gene_selector = widgets.SelectMultiple(options=self.var_names, description='Genes')
        self.cell_type_selector = widgets.SelectMultiple(options=self.adata.obs['cell type'].unique(), description='Cell Types')
        self.scale_widget = widgets.Dropdown(options=['var', 'obs', None], value='var', description='Scale:')
        self.color_map_widget = widgets.Dropdown(options=['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'Blues'], value='Blues', description='Color Map:')
        self.update_button = widgets.Button(description='Update Graph')
        self.output_widget = widgets.Output()

        # Configuration of events
        self.update_button.on_click(self.on_update_button_clicked)

    # Function for graph button 
    def on_update_button_clicked(self, b):
        selected_genes = self.gene_selector.value
        selected_cell_types = self.cell_type_selector.value
        scale = self.scale_widget.value
        color_map = self.color_map_widget.value

        with self.output_widget:
            self.output_widget.clear_output()
            self.var_names = list(selected_genes)
            self.selected_cell_type = list(selected_cell_types)
            self.show(scale=scale, cmap=color_map)

    # Function for graphic 
    def show(self, scale=None, cmap='Blues'):
        if self.var_names:
            filtered_data = self.adata[self.adata.obs['cell type'].isin(self.selected_cell_type)]
            sc.pl.matrixplot(filtered_data, self.var_names, 'cell type', 
                             dendrogram=True, cmap=cmap, standard_scale=scale, 
                             colorbar_title='column scaled\nexpression')

    # Function to display widgets
    def display_widgets(self):
        display(self.gene_selector, self.cell_type_selector, self.scale_widget, self.color_map_widget, self.update_button, self.output_widget)

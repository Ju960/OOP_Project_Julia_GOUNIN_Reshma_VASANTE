# pip install scanpy
# pip install igraph
# pip install leidenalg
import scanpy as sc
import ipywidgets as widgets
from IPython.display import display
from scanpy.plotting._matrixplot import MatrixPlot
import os
import matplotlib.pyplot as plt


class InteractiveMatrixPlot(MatrixPlot):
    def __init__(self, adata, var_names, **kwargs):
        # Filter genes to ensure they are present in adata
        filtered_var_names = [gene for gene in var_names if gene in adata.var_names]
        super().__init__(adata, filtered_var_names, groupby='cell type', **kwargs)
        self.adata = adata
        self.var_names = var_names

        # Creation of widgets
        self.gene_selector = widgets.SelectMultiple(options=filtered_var_names, description='Genes')
        self.cell_type_selector = widgets.SelectMultiple(options=adata.obs['cell type'].unique(), description='Cell Types')
        self.color_map_widget = widgets.Dropdown(options=['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'Blues'], value='Blues', description='Color Map:')
        self.update_button = widgets.Button(description='Update Graph')
        self.export_button = widgets.Button(description='Export Matrix as PNG')
        self.output_widget = widgets.Output()

        # Configuration of events
        self.update_button.on_click(self.on_update_button_clicked)
        self.export_button.on_click(self.on_export_button_clicked)

    # Function for graph button 
    def on_update_button_clicked(self, b):
        selected_genes = self.gene_selector.value
        selected_cell_types = self.cell_type_selector.value
        color_map = self.color_map_widget.value

        with self.output_widget:
            self.output_widget.clear_output()
            self.var_names = list(selected_genes)
            self.selected_cell_type = list(selected_cell_types)
            self.show_matrixplot(color_map=color_map)

    # Function for export button
    def on_export_button_clicked(self, b):
        export_filename = widgets.Text(value='matrix_plot.png', description='Filename:')
        export_button = widgets.Button(description='Export')
        export_output = widgets.Output()

        # Function for handling export button click
        def export_callback(button):
            with export_output:
                export_output.clear_output()
                filename = export_filename.value.strip()
                if not filename.endswith('.png'):
                    filename += '.png'

                # Check if selected genes are present in the Anndata object
                if all(gene in self.var_names for gene in self.gene_selector.value):
                    self.export_matrix_as_png(filename)
                    print(f'Matrix plot saved as {filename}')
                else:
                    print("Error: Selected genes not found in the Anndata object.")

        export_button.on_click(export_callback)

        # Display export widgets
        display(widgets.VBox([export_filename, export_button, export_output]))

    def export_matrix_as_png(self, filename='matrix_plot.png'):
        selected_genes = self.gene_selector.value
        selected_cell_types = self.cell_type_selector.value
        color_map = self.color_map_widget.value

        #Filter genes to ensure they are in adata.var_names
        valid_genes = [gene for gene in selected_genes if gene in self.adata.var_names]

        #Verify if valid genes are selected
        if valid_genes:
            #Create a directory to save images if it does not exist
            output_dir = 'output_images'
            os.makedirs(output_dir, exist_ok=True)

            #Construct the path of the output file
            output_path = os.path.join(output_dir, filename)

            #Filter the adata for the selected cell types
            filtered_data = self.adata[self.adata.obs['cell type'].isin(selected_cell_types)]
        
            #Create the plot with the valid genes
            sc.pl.matrixplot(filtered_data, valid_genes, 'cell type',
                         dendrogram=True, cmap=color_map,
                         colorbar_title='column scaled\nexpression', show=False)

            #Save the figure
            plt.savefig(output_path, bbox_inches='tight', dpi=300)
            plt.close()
            print(f'Graphique enregistré sous : {output_path}')
        else:
            #Display an error message if no valid genes are selected
            print("Erreur : Aucun gène valide sélectionné pour l'exportation.")


    #Function to display widgets
    def display_widgets(self):
        display(self.gene_selector, self.cell_type_selector, self.color_map_widget,
                self.update_button, self.export_button, self.output_widget)

    def show_matrixplot(self, color_map='Blues'):
        valid_genes = [gene for gene in self.var_names if gene in self.adata.var_names]
    
        if valid_genes:
            filtered_data = self.adata[self.adata.obs['cell type'].isin(self.selected_cell_type)]

            sc.pl.matrixplot(filtered_data, valid_genes, 'cell type',dendrogram=True, cmap=color_map, 
                         colorbar_title='column scaled\nexpression')
        else:
            print("Aucun gène valide sélectionné.")

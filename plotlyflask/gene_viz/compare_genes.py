
import numpy as np

from plotlyflask.gene_viz import shared
import plotly.graph_objects as go

def plot(found_in, metadata, selected_genes, datasets_selected, plot_type, xaxis_type):
    if not found_in.empty:
        ret_str = 'Genes {} were found'.format(selected_genes)

        gene_A, gene_B = found_in["genes"].values[0], found_in["genes"].values[1]

        # create the figure layout
        layout = go.Layout(
            title=" {} vs {} ".format(gene_A, gene_B),
            height=700
        )
        metadata_selected = shared.filter_data(metadata, "Dataset", datasets_selected)

        prcsd_df = found_in.iloc[0]
        prcsd_df = prcsd_df.drop(["genes"]).T.reset_index()
        prcsd_df.columns = ["Sample", "TPM"]
        prcsd_df = shared.sync_df(metadata_selected, prcsd_df)
        prcsd_df.rename(columns={"TPM": gene_A}, inplace=True)
        prcsd_df[gene_B] = found_in[prcsd_df["Sample"].values].iloc[1].values

        hover_data = ["Sample", "Tissue", "Dataset", "subset_name",
                      "NHU_differentiation", "Gender", "TER", "Substrate"]

        config = {"x": gene_A, "y": gene_B, "color": "Dataset",
                  "hover_data": hover_data, "xaxis_type": xaxis_type}
        fig = shared.select_plot_type(prcsd_df, config, plot_type, isComparison=True)

        fig.update_layout(layout)

        if xaxis_type != "log10":
            offset_x = prcsd_df[gene_A].values.max() * 0.015
            offset_y = prcsd_df[gene_B].values.max() * 0.015
            fig.update_xaxes(
                range=[-offset_x, prcsd_df[gene_A].values.max() + offset_x])
            fig.update_yaxes(
                range=[-offset_y, prcsd_df[gene_B].values.max() + offset_y])
        else:
            offset_x = np.log10(list(prcsd_df[gene_A].values)).max() * 0.015
            offset_y = np.log10(list(prcsd_df[gene_B].values)).max() * 0.015
            fig.update_xaxes(
                range=[-offset_x, np.log10(list(prcsd_df[gene_A].values)).max() + offset_x])
            fig.update_yaxes(
                range=[-offset_y, np.log10(list(prcsd_df[gene_B].values)).max() + offset_y])

        return ret_str, fig, prcsd_df

    return "Genes {} were not found in any of the datasets".format(selected_genes), {"data": {}}, found_in

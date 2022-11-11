import plotly.graph_objects as go

from plotlyflask.gene_viz import shared


def plot(found_in, metadata, user_input, datasets_selected, ter, plot_type, xaxis, xaxis_type):
    if not found_in.empty:
        ret_str = 'Gene {} was found'.format(user_input)

        # create the figure layout
        layout = go.Layout(
            title="Swarm plot for {}".format(found_in["genes"].values[0]),
            height=700
        )
        # prepare datafrarme
        found_in = found_in.drop(["genes"], axis=1).T.reset_index()
        found_in.columns = ["Sample", "TPM"]

        # process metadata and add to the df
        metadata_selected = shared.filter_data(metadata, "Dataset", datasets_selected)
        processed_df = shared.sync_df(metadata_selected, found_in)

        hover_data = ["Sample",  "Dataset", "subset_name", "shared_num_same_col",
                      "Tissue", "NHU_differentiation", "Gender", "TER", "Substrate"]
        if not processed_df.empty:
            if not datasets_selected:
                return "Select a dataset first", {"data": {}}, found_in

            config = {"x": xaxis, "y": "TPM", "color": "Dataset",
                      "hover_data": hover_data, "xaxis_type": xaxis_type, }

            # create the main trace
            fig = shared.select_plot_type(processed_df, config,
                                   plot_type, isComparison=False, ter=ter)

            if xaxis_type == "linear":
                if processed_df["TPM"].values.max() <= 25:
                    fig.update_yaxes(range=[0, 25])
                else:
                    offset_y = processed_df["TPM"].values.max() * 0.015
                    fig.update_yaxes(
                        range=[-offset_y, processed_df["TPM"].values.max() + offset_y])

            fig.update_layout(layout)
            return ret_str, fig, processed_df

    return "Gene {} not found in any of the datasets".format(user_input), {"data": {}}, found_in

# Others


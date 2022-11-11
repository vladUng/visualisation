
import plotly.express as px


def select_plot_type(df, config, plot_type, isComparison=False, ter=None):
    x, y, color, hover_data, xaxis_type = config["x"], config[
        "y"], config["color"], config["hover_data"], config["xaxis_type"]
    ret_fig = {}
    categorical_order = []

    log = False
    if xaxis_type != "linear":
        log = True
        if not isComparison:
            df.loc[df[df["TPM"] < 1].index.values]["TPM"] = 1.0
        else:
            df.loc[df[df[y] < 1].index.values][y] = 1.0
            df.loc[df[df[x] < 1].index.values][x] = 1.0

    if not isComparison:
        # slightly modified the subtypes and duplicate the ones that have tight TER barrier
        if ter != None and "tight" in ter:
            tight_limit = 500
            ter_col = []
            for _, row in df.iterrows():
                new_name = ""
                if row["isTER"]:
                    if float(row["TER"]) >= tight_limit:
                        new_name = row[x]+"_tight"
                    else:
                        new_name = row[x]+"_non-tight"
                else:
                    new_name = row[x]

                ter_col.append(new_name)

                if new_name not in categorical_order:
                    categorical_order.extend([new_name])

            categorical_order.sort()
            df[x] = ter_col
            color = x

        # do the plot
        if plot_type == "violin":
            ret_fig = px.violin(df, x=x, y=y, color=color,
                                box=True, hover_data=hover_data, log_y=log)
        elif plot_type == "violin_points":
            ret_fig = px.violin(df, x=x, y=y, color=color, box=True,
                                points="all", hover_data=hover_data, log_y=log)
        elif plot_type == "box":
            ret_fig = px.box(df, x=x, y=y, color=color, hover_data=hover_data,
                             log_y=log, category_orders={x: categorical_order})
            ret_fig.update_traces(boxmean="sd")
        elif plot_type == "box_points":
            ret_fig = px.box(df, x=x, y=y, color=color, points="all",
                             hover_data=hover_data, log_y=log,  category_orders={x: categorical_order})
            ret_fig.update_traces(boxmean="sd")
        else:
            ret_fig = px.strip(df, x=x, y=y, color=color, hover_data=hover_data,
                               log_y=log, category_orders={x: categorical_order})
    else:
        ret_fig = px.strip(df, x=x, y=y, color=color,
                           hover_data=hover_data, log_y=log, log_x=log)

    return ret_fig

def sync_df(df_metadata, df):
    """ The role of this functions is to make sure that metadata df and all_data corresponds to the same samples after we've filter out some off the samples.

    Args:
        df_metadata ([DataFrame]): The medata
        df ([DataFrame]): TPM dataframe

    Returns:
        DataFrame: The metadata sorted
    """
    # the df with common samples from both df metadata and the df with tpms
    df_common = filter_data(df, "Sample", list(df_metadata["Sample"].values))

    # for cases when some samples are in a database but not others
    metadata_common = filter_data(df_metadata, "Sample", list(df["Sample"].values))

    # sort both by samples so that there are syncd
    metadata_sorted = metadata_common.sort_values("Sample")
    df_common = df_common.sort_values("Sample")
    metadata_sorted["TPM"] = df_common["TPM"].values
    return metadata_sorted

def filter_data(df, col, values_to_keep):
    if values_to_keep:
        df = df[df[col].isin(values_to_keep)]
    return df

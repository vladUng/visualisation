import numpy as np
import pandas as pd

from plotlyflask.gene_diff import gene_markers as gm 

import time

from os import path

def create_custom_traces(selected_genes = None):
    custom_traces = []
    if selected_genes:
        custom_traces.append({"genes":selected_genes, "title": "Selected genes"})

    ifnq_genes = ["INFG", 'B2M', 'BATF2', 'BTN3A3', 'CD74', 'CIITA', 'CXCL10', 'CXCL9','EPSTI1', 'GBP1', 'GBP4', 'GBP5', 'HAPLN3', 'HLA-DPA1', 'IDO1', 'IFI30', 'IFI35', 'IFIT3', 'IFITM1', 'IL32', 'IRF1', 'NLRC5', 'PARP9', 'PSMB8-AS1', 'PSMB9', 'SAMHD1', 'SECTM1', 'STAT1', 'TAP1','TNFSF13B', 'TYMP', 'UBE2L6', 'WARS1', 'ZBP1']

    diff_neuronal = ['ST18', 'DPYSL5', 'DCX', 'TMEM145', 'ELAVL3', 'TOP2A', 'SCN3A', 'ASXL3', 'SEZ6', 'ATP1A3', 'CELF3', 'PEX5L', 'PCSK2', 'PROX1', 'EPHA7', 'CHGB', 'UNC13A', 'RIMS2', 'CAMK2B', 'NKX2-2', 'AP3B2']

    # diff when Consensus vs Mixed and LumInf vs Mixed
    ne_dif = ["ARHGAP33","CBX5","CCNE2","CDCA5","CDCA8","CDK5R1","CDKAL1","CENPF","CLSPN","CNIH2","COCH","DCX","DPYSL5","E2F7","ENHO","FBXO43","FBXO5","FSD1","GNAZ","GPRIN1","GTSE1","HES6","HMGB2","IFT81","KCNH3","LHX2","LMNB1","MAD2L1","MCM2","MSANTD3-TMEFF1","MT3","PBK","PDXP","PIMREG","PPM1D","PRAME","PSIP1","PSRC1","PTTG1","QSOX2","RCOR2","SALL2","SCML2","SMC2","SRSF12","STMN1","TEDC2","TLCD3B","TMSB15A","TUBA1B","TUBB2B","TXNDC16","TYMS","USP1","WASF1","ZNF711", "C12orf75","CBX2","CDC25C","DEPDC1","EZH2","GAS2L3","KIF18B","NUSAP1","ODC1"]

    mixed_lumInf_diff =["ADGRF4","ADIRF","ALS2CL","ANXA8","ANXA8L1","B3GNT3","BATF","BCAS1","C10orf99","C19orf33","CH17-360D5.2","CLDN1","CTSH","CXCL17","EPS8L1","EVPL","FAM110C","FXYD3","GBP2","GNA15","GPR87","IL1RN","ITGB4","ITGB6","KCNN4","KRT19","KRT7","MPZL2","MYOF","P2RY2","PLEKHN1","PRSS22","PSCA","PTGES","PTK6","S100A11","S100A6","S100P","SDC1","SNCG","SYT8","SYTL1","TACSTD2","TINAGL1","TMEM40","UPK1A","UPK2","UPK3A","UPK3B","VGLL1","WNT7B", "ACSF2","ARHGAP27","BICDL2","CAPS","CARD11","CBLC","CLDN4","CSTB","CYP4B1","DENND2D","DTX4","EHF","ELF3","EPHA1","EPN3","FBP1","FOXQ1","GATA3","GDF15","GGT6","GPR160","GRHL1","GRHL3","IQANK1","KRT7-AS","LLNLR-231D4.1","LPAR5","METRNL","NECTIN4","OVOL1","PLA2G2F","PLA2G4F","PPL","PROM2","PSD4","RAB25","RBBP8NL","RBM47","RP4-568C11.4","S100A14","SCNN1B","SEMA4A","SPINK1","SSH3","TFAP2C","TJP3","TMC6","TMPRSS2","TRPM4","UGT1A6","VAMP8","VSIG2"]

    jens_work = ["GJB1"]

    basal_small_specific = ['AC005301.9',   'AIF1', 'APOC2', 'CALHM6', 'CCL4', 'CD163',  'CLCF1',  'DERL3',  'ELN',  'FBLN2',  'FCGR2A',  'FOLR2',  'FST',  'GZMB',  'HCST',  'IFIT3',  'ITGB2',  'KRT80',  'MFAP4',  'MMP23B',  'MS4A6A',  'MT1M',  'MZB1',  'NKG7',  'NR4A1AS',  'OLFML2A',  'OLFML3',  'RGS2',  'RP11-54H7.4',  'TRAC']

    custom_traces.append({"genes": basal_small_specific, "title": "Small Ba/Sq specific"})

    # custom_traces.append({"genes": jens_work, "title": "Jen's work"})

    # custom_traces.append({"genes":ne_dif, "title": "Diff for NE"})

    # custom_traces.append({"genes":mixed_lumInf_diff, "title": "Diff for Mixed/LumInf"})


    # SB
    # custom_traces.append({"genes":ifnq_genes, "title": "SB_IFNQ"})
    # custom_traces.append({"genes":diff_neuronal, "title": "Diff old vs remap"})

    # ryan_genes = ["FGFR3", "EGFR", "TP53"]
    custom_traces = gm.add_tcga_markers(custom_traces=custom_traces)
    # custom_traces = gm.add_lund_markers(custom_traces=custom_traces)

    # custom_traces = gm.lumInf_mixed(custom_traces=custom_traces)

    custom_traces = gm.significant_genes(custom_traces=custom_traces)

    custom_traces = gm.low_significant_genes(custom_traces=custom_traces)
    
    return custom_traces
    # return []

def create_gene_trace(df, genes, name="custom genes", marker_color="yellow", marker_size=12, df_2=None): 

    selected_df = df[df["genes"].isin(genes)]

    x, y = [], []
    if df_2 is None:
        y = -np.log10(selected_df["q"])
        x = selected_df['fold_change']
    else:
        selected_df_2 = df_2[df_2["genes"].isin(genes)]
        
        x = -np.log10(selected_df["q"]) * selected_df["fold_change"]
        y = -np.log10(selected_df_2["q"]) * selected_df_2["fold_change"]

    markers = {"size": marker_size, "color": selected_df.shape[0] * [marker_color], "symbol": "x" }

    trace = dict(type='scatter', x=x, y= y,  showlegend=True, marker=markers, text=selected_df["genes"], mode="markers+text" , name=name, textposition= 'top center')

    return trace

##### Functions for .tsv handling #####
def import_data(fullPath):
    start_time = time.time()
    if path.exists(fullPath):
        ret_dict = {}
        ret_dict["data"] = pd.read_csv(fullPath, delimiter="\t")
        print("Finished loading the data in {}".format(time.time()-start_time))
        return ret_dict
    else:
        return None

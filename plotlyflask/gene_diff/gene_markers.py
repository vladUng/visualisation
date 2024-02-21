#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   gene_markers.py
@Time    :   2022/10/27 13:18:54
@Author  :   Vlad Ungureanu
@Version :   1.0
@Contact :   vlad.ungureanu@york.ac.uk
@Desc    :   File to put genes that are show in the DEA/Pi/Scatter plot
'''

def add_uroth_bladder_cancer(custom_traces):
    """

    These genes were selected by Ryan and are out of interests in bladder/urothelium 

    Args:
        custom_traces (_type_): _description_
    """

    tf_diff = ["P63", "FOXA1","PPARG", "RARG", "IRF1", "ELF3", "GRHL3", "KLF5", "GATA4", "GATA6", "GATA3"]
    egfr_fam = ["EGFR", "ERBB2", "ERBB3", "ERBB4", "EGF", "AREG", "HBEGF","TGFA","BTC", "EREG"]
    fgfr_fam = ["FGFR1", "FGFR2", "FGFR3", "FGF1", "FGF2"]
    map_kpathway = ["RAS", "RAF", "MEK1", "MEK2", "MEK3", "MEK4","ERK"]
    pi3_kpathway = ["PIK3C3", "PIK3R2", "PIK3C2B", "AKT1", "AKT2"]
    krt = ["KRT13", "KRT14", "KRT15", "KRT20"]
    cld = ["CLDN3", "CLDN4", "CLDN5" ]
    others = ["MKI67", "MCM2", "UPK3A", "ZO1", "TJP1", "ZO2", "TJP2", "ZO3", "TJP3"]
    hox_ur = ["HOXB2", "HOXB3", "HOXB5", "HOXB6", "HOXB8"]
    hox_bla = ["HOXA9", "HOXA10", "HOXA11", "HOXA13"]


    custom_traces.append({"genes": tf_diff, "title": "TF Diff"})
    custom_traces.append({"genes": egfr_fam, "title": "EGFR fam"})
    custom_traces.append({"genes": fgfr_fam, "title": "FGFR fam"})
    custom_traces.append({"genes": map_kpathway, "title": "MAP kin"})
    custom_traces.append({"genes": pi3_kpathway, "title": "PI3 kin"})
    custom_traces.append({"genes": krt, "title": "Diff markers"})
    custom_traces.append({"genes": cld, "title": "Urothelial diff"})
    custom_traces.append({"genes": hox_ur, "title": "Hox ureteric"})
    custom_traces.append({"genes": hox_bla, "title": "Hox bladder"})

    return custom_traces


def gc42_high_genes(custom_traces):

    lump = ['AQP3', 'IGLC2', 'DHRS2', 'MAOA', 'TMPRSS4', 'GPX2', 'HPGD', 'VSIG2', 'SPINK1', 'TBX3', 'PSCA', 'UPK2']
    lumInf_ns = ['MUC20', 'UPK2', 'H19', 'CTSK', 'VGLL1', 'IGKV3-20', 'PSCA', 'IGLC2', 'IGHJ4', 'AEBP1', 'KRT23', 'GPX2', 'IGLC1', 'POSTN', 'IGHG1', 'ACTG2', 'IGHV3-23', 'IGHA1', 'IGLC3']
    large_basal = ['AQP3', 'TNS4', 'KRT16', 'MALL', 'LY6D', 'SLPI', 'FAM83A', 'DSP', 'KRT5', 'SERPINB5', 'LCN2', 'S100A8', 'COL17A1', 'CDH3', 'ANXA8L1']
    small_basal =['IGKJ2', 'IGKV1-5', 'HLA-DQA1', 'IGKV1D-39', 'CTHRC1', 'H19', 'KRT5', 'CTSK', 'FCGR3A', 'IGKV3-15', 'IGKJ1', 'IGHJ5', 'CDH3', 'AIF1', 'IGKV3-20', 'IGHG3', 'IGLC2', 'IGHJ4', 'AEBP1', 'IGKV4-1', 'DSP', 'IGKV2-28', 'IGLC1', 'IGLV2-14', 'POSTN', 'IGHG1', 'IGHG2', 'IGHV3-23', 'IGLC3', 'IGKJ4', 'LYZ']
    ne = ['H2AC19', 'IGLC2', 'IGHG1', 'IGFBP2', 'H19']

    custom_traces.append({"genes": lump, "title": "High LumP"})
    custom_traces.append({"genes": lumInf_ns, "title": "High LumInf/NS"})
    custom_traces.append({"genes": large_basal, "title": "High Large Ba/Sq"})
    custom_traces.append({"genes": small_basal, "title": "High Small Ba/Sq"})
    custom_traces.append({"genes": ne, "title": "High NE"})

    return custom_traces

def gc42_low_genes(custom_traces):

    lump = ['ENSG00000289474', 'ENSG00000291105', 'ASPN', 'LMOD1', 'SUGCT', 'TGM1', 'CES1', 'ITM2A', 'IL20RB', 'S100A3', 'ANPEP', 'MIR663AHG', 'MT1M', 'PTHLH', 'TSPAN7', 'APLP1', 'SPP1', 'SYNM', 'GZMA', 'SELL', 'HSD17B1', 'GZMB', 'MEF2C', 'HSDL2-AS1', 'EFHD1', 'SOX2-OT', 'C4A', 'WNT5B', 'F13A1', 'IGHM', 'MMP9', 'SCARNA9']
    lumInf_ns =  ['ENSG00000269968', 'GPC3', 'TGM1', 'CES1', 'ID4', 'ITM2A', 'FABP3', 'IL20RB', 'S100A3', 'ENSG00000251034', 'ENSG00000286746', 'ANPEP', 'MIR663AHG', 'ENSG00000274281', 'CFLAR-AS1', 'PTHLH', 'ENSG00000289609', 'ENSG00000267458', 'GZMA', 'ENSG00000290560', 'LINC01088', 'SCARNA9', 'CRMP1', 'ENSG00000286125', 'ENSG00000237356', 'ENSG00000283633', 'BCL10-AS1', 'HSDL2-AS1', 'ALDH1A1', 'SCGB3A2', 'ENSG00000232035', 'SOX2-OT', 'ENSG00000275216', 'MSX1', 'ENSG00000286827', 'HES6', 'CD36', 'LINC03000', 'RASD1', 'MMP9', 'ENSG00000291060', 'ENSG00000230612']
    large_basal =['H2BC18', 'ENSG00000269968', 'LMOD1', 'ID4', 'CES1', 'FABP3', 'ITM2A', 'NRN1', 'SORBS1', 'ENSG00000288934', 'PLA2G2A', 'S100A3', 'ENSG00000251034', 'ENSG00000286746', 'ENSG00000274281', 'HSPB6', 'APLP1', 'TSPAN7', 'ENSG00000267458', 'SYNM', 'FXYD1', 'SELL', 'ENSG00000287516', 'AIFM3', 'HSDL2-AS1', 'EFHD1', 'MYH11', 'ALDH1A1', 'ENSG00000232035', 'TESC', 'SOX2-OT', 'ENSG00000286827', 'HES6', 'CD36', 'F13A1', 'RASD1', 'MMP9', 'SCARNA9']
    small_basal = ['H2BC18', 'ENSG00000269968', 'LINC02895', 'GPC3', 'TGM1', 'CES1', 'FABP3', 'FCRLB', 'NRN1', 'GPAM', 'IL20RB', 'S100A3', 'ENSG00000251034', 'ENSG00000286746', 'AKR1C1', 'ZBTB20', 'ENSG00000274281', 'MIR663AHG', 'CFLAR-AS1', 'ENSG00000289609', 'SPINK1', 'ENSG00000267458', 'ENSG00000291060', 'ENSG00000290560', 'SELL', 'LINC01088', 'HSD17B1', 'ENSG00000237356', 'ENSG00000286125', 'ENSG00000287516', 'AIFM3', 'HSDL2-AS1', 'ALDH1A1', 'ENSG00000232035', 'TESC', 'SOX2-OT', 'MSX1', 'SPRR1A', 'ENSG00000286827', 'HES6', 'CD36', 'ENSG00000291168', 'DANT2', 'FABP4', 'RASD1', 'CYP3A5', 'SCIN', 'SCARNA9', 'GSTM2', 'ENSG00000230612']
    ne = ['ENSG00000269968', 'ENSG00000289474', 'LMOD1', 'SUGCT', 'TGM1', 'FABP3', 'ITM2A', 'PLA2G2A', 'OLR1', 'RRAD', 'ANPEP', 'MIR663AHG', 'CFLAR-AS1', 'ENSG00000267458', 'SYNM', 'TSPAN7', 'GZMA', 'KLHDC7B', 'CCL20', 'LINC01088', 'SCARNA9', 'S100A8', 'GZMB', 'ENSG00000287516', 'MYH11', 'CST6', 'PDZK1IP1', 'LASTR', 'TNNI2', 'ENSG00000232035', 'ENSG00000291060']

    custom_traces.append({"genes": lump, "title": "Low LumP"})
    custom_traces.append({"genes": lumInf_ns, "title": "Low LumInf/NS"})
    custom_traces.append({"genes": large_basal, "title": "Low Large Ba/Sq"})
    custom_traces.append({"genes": small_basal, "title": "Low Small Ba/Sq"})
    custom_traces.append({"genes": ne, "title": "Low NE"})

    return custom_traces

def gc42_uniq_genes(custom_traces):

    small_basal=['ENSG00000291105', 'OLFML2A', 'APOC2', 'ANGPTL2', 'CALHM6', 'GLIPR2', 'CD163', 'TRBC1', 'ANPEP', 'FCGR3A', 'FCGR1A', 'OLFML3', 'PDGFRB', 'VCAN', 'FCGR2A', 'ITGB2', 'ENSG00000289273', 'TRAC', 'GZMB', 'HCST', 'NKG7', 'CCL4', 'CCL5', 'MS4A6A', 'VSIG4']
    ne=['IGSF9', 'SNCA', 'DMTN', 'PRKDC', 'ZNF331', 'TDRP', 'DSC2', 'MXRA8', 'ZSCAN18', 'ADAMTSL2', 'DNM1', 'ARHGAP4', 'FGFRL1', 'SOX9', 'RCAN2', 'HIC2', 'TSPAN33', 'KITLG', 'IL17D', 'MLLT11', 'ST6GAL1', 'KHDRBS3', 'CADM1', 'RAP1GAP', 'LIPE']
    lump=['PYROXD2', 'HOXB4', 'ENSG00000272141', 'ENSG00000269968', 'FOXA1', 'GPC3', 'KLRK1-AS1', 'DUSP2', 'SPIRE2', 'BHLHE41', 'FRY', 'ATF7IP2', 'HSH2D', 'EPG5', 'SLC40A1', 'ALDH1A2', 'HSD17B14', 'TBC1D3L', 'TMEM184A', 'CEACAM1', 'ABCC3', 'HOXB5', 'CXCL12', 'SPOCD1', 'HOXB2']
    large_basal=['FKBP5', 'CDK5R1', 'ANO1', 'KLC3', 'ENSG00000290574', 'LRIG3', 'TGM1', 'MAFB', 'CSPG4', 'CLEC2B', 'TRIM22', 'S1PR5', 'IL20RB', 'HLA-DOA', 'OLR1', 'FAP', 'IL15RA', 'GJB5', 'ADAM19', 'LINC01013', 'TRIM7', 'DCBLD2', 'STC1', 'SLC38A5', 'HS3ST1']
    lumInf_ns=['QPCT', 'EDARADD', 'PTGDS', 'ASPN', 'BEX2', 'C3orf14', 'PRSS27', 'CAPN12', 'ENTPD3', 'ALDH3B1', 'ANXA9', 'SUSD2', 'GATM', 'FHL1', 'DAB2', 'SECTM1', 'RDX', 'TMEM125', 'SELL', 'ACER2', 'FADS1', 'RHPN2', 'SPOCK2', 'CFD', 'NCF4']

    custom_traces.append({"genes": lump, "title": "Uniq LumP"})
    custom_traces.append({"genes": lumInf_ns, "title": "Uniq LumInf/NS"})
    custom_traces.append({"genes": large_basal, "title": "Uniq Large Ba/Sq"})
    custom_traces.append({"genes": small_basal, "title": "Uniq Small Ba/Sq"})
    custom_traces.append({"genes": ne, "title": "NE uniq"})

    return custom_traces 

def add_tcga_markers(custom_traces):

    luminal_markers = ["KRT20", "PPARG", "FOXA1", "GATA3", "SNX31", "UPK1A", "UPK2", "FGFR3"]
    basal_markers = ["CD44", "KRT6A", "KRT5", "KRT14", "COL17A1"]
    squamos_markers = ["DSC3", "GSDMC", "TCGM1", "PI3", "TP63"]
    immune_markers = ["CD274", "PDCD1LG2", "IDO1", "CXCL11", "L1CAM", "SAA1"]
    neural_diff = ["MSI1", "PLEKHG4B", "GNG4", "PEG10", "RND2", "APLP1", "SOX2", "TUBB2B"]

    # TCGA
    custom_traces.append({"genes":luminal_markers, "title": "TCGA_luminal"})
    custom_traces.append({"genes":basal_markers, "title": "TCGA_basal"})
    custom_traces.append({"genes":immune_markers, "title": "TCGA_immune"})
    custom_traces.append({"genes":squamos_markers, "title": "TCGA_squamos"})
    custom_traces.append({"genes":neural_diff, "title": "TCGA_neuroendocrine"})

    return custom_traces

# lund
def add_lund_markers(custom_traces):
    # both ba_sq_inf and mes_like
    lund_qtc1 = ["FLI1", "FOXP3", "ILKZF1", "IRF4", "IRF8", "RUNX3", "SCML4", "SPI1", "STAT4", "TBX21", "TFEC"]
    lund_qtc2 = ["AEBP1", "BNC2", "GLI2", "GLIS1", "HIC1", "MSC", "PPRX1", "PPRX2", "TGFB1I1", "TWIST1"]
    lund_qtc3 = ["EBF1", "HEYL", "LEF1", "MEF2C", "TCF4", "ZEB1", "ZEB2"]
    lund_qtc8 = ["GATA5", "HAND1", "HAND2", "KLF16"]
    lund_qtc17 = ["ARID5A", "BATF3", "VENTX"]
    lund_ba_mes = lund_qtc1 + lund_qtc2 + lund_qtc3 + lund_qtc8 + lund_qtc17

    lund_ba_sq = ["BRIP1", "E2F7", "FOXM1", "ZNF367", "IRF1", "SP110", "STAT1"]
    lund_mes = ["TP53", "RB1", "FGFR3", "ANKHD1", "VIM", "ZEB2"]
    ba_sq_inf = ["CDH3", "EGFR"]

    lund_sc_ne = ["CHGA", "SYP", "ENO2", "EPCAM"] #Highly expressed

    custom_traces.append({"genes":lund_ba_mes, "title": "Mes-like + Ba/sq-inf"})
    custom_traces.append({"genes":lund_ba_sq, "title": "Lund ba/sq"})
    custom_traces.append({"genes":lund_mes, "title": "Mes-like"})
    custom_traces.append({"genes":ba_sq_inf, "title": "Ba/sq-inf"})

    custom_traces.append({"genes": lund_sc_ne, "title": "Sc/NE"})
    
    return custom_traces

# immune cells
def add_imune_markers(custom_traces):
    b_cells = ["BCL2", "BCL6", "CD19", "CD1D", "CD22", "CD24", "CD27", "CD274","CD34", "CD38", "CD40","CD44","CD5","CD53","CD69","CD72", "CD79A", "CD79B", "CD80", "CD86", "CD93", "CR2", "CXCR4", 'CXCR5',"FAS","FCER2", "FCRL4" "HAVCR1","IL10", 'IL2RA','IL7R','IRF4','ITGAX', 'LILRB1','MME','MS4A1','NT5E','PDCD1LG2','PRDM1','PTPRC','SDC1','SPN','TFRC','TLR9','TNFRSF13B','TNFRSF13C','TNFRSF17','XBP1']
    t_cells = ['CD4', 'CD8', 'CCR4', 'CCR5', 'CCR6', 'CCR7', 'CCR10', 'CD127', 'CD27', 'CD28', 'CD38', 'CD58', 'CD69', 'CTLA4', 'CXCR3', 'FAS', 'IL2RA',
            'IL2RB', 'ITGAE', 'ITGAL', 'KLRB1', 'NCAM1', 'PECAM1', 'PTGDR2', 'SELL', 'IFNG', 'IL10', 'IL13', 'IL17A', 'IL2', 'IL21','IL22', 'IL25', 'IL26', 'IL4', 'IL5', 'IL9', 'TGFB1', 'TNF', 'AHR', 'EOMES','FOXO4', 'FOXP1', 'FOXP3', 'GATA3','IRF4', 'LEF1', 'PRDM1', 'RORC','STAT4', 'TBX21','TCF7', 'GZMA']

    nk_cells = ['B3GAT1','CCR7','CD16','CD2','CD226','CD244','CD27','CD300A','CD34','CD58','CD59','CD69','CSF2','CX3CR1','CXCR1','CXCR3','CXCR4','EOMES','GZMB','ICAM1','IFNG','IL1R1','IL22','IL2RB','IL7R','ITGA1','Itga2','ITGAL','ITGAM','ITGB2','KIR2DL1','KIR2DL2','KIT','Klrb1c','KLRC1','KLRC2','KLRD1','KLRF1','KLRG1','KLRK1','LILRB1','Klra4','Klra8','NCAM1','NCR1','NCR2','NCR3','PRF1','SELL','SIGLEC7','SLAMF6','SPN','TBX21','TNF']

    macrophages_cells = [ 'ADGRE1','CCR2','CD14','CD68','CSF1R','Ly6c1','MARCO','MRC1','NOS2','PPARG','SIGLEC1','TLR2','ARG1','CD163','CD200R1','CD80','CD86','CLEC10A','CLEC7A','CSF2','CX3CR1','FCGR1A','ITGAM','MERTK','PDCD1LG2','Retnla','TNF','CCL22','CD36','CD40','IL10','IL1B','IL6','LGALS3','TLR4','CCL2','CCR5','CD209','CD63','CD86','CSF1','CXCL2','FCGR3A','IFNG','IL4','IRF4','ITGAX','MSR1','PDGFB','PTPRC','STAT6','TIMD4','Chil3','CLEC6A','IL1R1','ITGB2','PDCD1LG2','TLR7']

    monocyte_cells = ['CD14','CD16','CSF1R','CX3CR1','ITGAM','ITGAX','LY6C1','CCR2','CXCR4','FCGR1A','SELL','SPN','ADGRE1','CCR7','TNF','CD86','IL10','IL1B','MERTK','TREML4','CD209','NR4A1','Ly6a','PTPRC','IL3RA','CD27','CCR5','CD32','CD1A','MRC1','ITGB3','CD9','CXCR6','CCR1','FLT3','KLF2','CLEC12A','CCR6','CCR8','CD68','CLEC7A','KIT','MAF','MAFB','SPI1','CD1C','PPARG','CEBPB','ITGAE','TEK']

    custom_traces.append({"genes":b_cells, "title": "B markers"})
    custom_traces.append({"genes":t_cells, "title": "T markers"})
    custom_traces.append({"genes":nk_cells, "title": "NK markers"})
    custom_traces.append({"genes":macrophages_cells, "title": "Macrophae markers"})
    custom_traces.append({"genes":monocyte_cells, "title": "Monocitye markers"})

    return custom_traces

def add_ig(custom_traces):
    interferon_g = ["IGHG1","IGHG3","IGHJ2","IGHJ3","IGHJ4","IGHM","IGHV1-69D","IGHV2-70","IGHV3-11","IGHV3-15","IGHV4-39","IGKC","IGKJ1","IGKJ4","IGKV1-5","IGKV1-6","IGKV1-9","IGKV1D-39","IGKV2-28","IGKV4-1","IGLC1","IGLC2","IGLC3","IGLJ2","IGLJ3","IGLV1-40","IGLV1-44","IGLV3-19"]

    sig_gens_int_h = ["CMKLR1","CYTH4","FAM78A","GIMAP1","GIMAP5","IFFO1","IL10RA","LILRB4","NRROS","PLEKHO1","RASSF4","SLAMF8","SLC7A7","TCN2","TFEC"]

    custom_traces.append({"genes":interferon_g, "title": "IG family "})
    custom_traces.append({"genes":sig_gens_int_h, "title": "Significant in Basal 2"})

    return custom_traces

def add_basal_specific(custom_traces):

    serp_genes = ['SERPINB13','SERPINB3', 'SERPINB4']
    spr_genes = ['SPRR1A', 'SPRR1B', 'SPRR2A', 'SPRR2D', 'SPRR3']

    # cell diff
    s100_genes = ['S100A2', 'S100A7', 'S100A8', 'S100A9']
    cell_diff = ["FGFBP1", 'HSPB1']
    keratines = ['KRT14', 'KRT17', 'KRT5','KRT6C']

    custom_traces.append({"genes":serp_genes, "title": "SERP family"})
    custom_traces.append({"genes":spr_genes, "title": "SPR family"})

    custom_traces.append({"genes":s100_genes, "title": "S100 family"})
    custom_traces.append({"genes":cell_diff + keratines, "title": "Cell diff + keratines"})

    high_exprs_genes_in_large_basal_vlad = ["AKR1C3","AQP3","ELF3","FABP5","GJB2","GJB5","GJB6","KLF5","KRT15","KRT16","KRT6B","LCN2","LYPD3","NDUFA4L2","NECTIN1","NECTIN4","PKP1","PROM2","PTHLH","RAB38","RHOV","S100A14","SEMA4B","SERPINB5","SLPI","THBD","TP63"]
    custom_traces.append({"genes":high_exprs_genes_in_large_basal_vlad, "title": "High expressed in basal large"})

    return custom_traces

# LumInf vs Mixed 
def lumInf_mixed(custom_traces):

    #highly expressed in LumInf
    lumInf_top_branch = ["AEBP1","C1R","C1S","C1orf162","C3","CCN2","COL6A3","COX7A1","CPA3","CTSK","DCN","DPT","EVI2A","FGF7","GAS7","GGT5","GLT8D2","GPR68","GREM1","HSPB1","IGFBP4","IRAG1","ISLR","LUM","LY96","MFAP4","MMP2","MRGPRF","NNMT","PDPN","PODN","RARRES2","RNASE6","SELPLG","SERPING1","TAGLN","THBS1","TIMP3","TMEM119","TMEM273","TPSAB1","TPSB2"]

    lumInf_middle_branch = ["ACTA2","ACTG2","ANXA8","ANXA8L1","C19orf33","C4A","C4B","CCL19","CCL21","CLIC3","CLU","CNN1","COL10A1","COL3A1","COMP","DES","HSPB6","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4","IGHJ2","IGHJ3","IGHJ4","IGHJ5","IGHM","IGHV1-2","IGHV1-3","IGHV1-46","IGHV1-69D","IGHV3-11","IGHV3-15","IGHV3-21","IGHV3-23","IGHV3-30","IGHV3-33","IGHV3-43","IGHV3-48","IGHV3-49","IGHV3-74","IGHV4-34","IGHV4-39","IGHV4-61","IGHV5-51","IGKC","IGKJ1","IGKJ2","IGKJ5","IGKV1-12","IGKV1-33","IGKV1-5","IGKV1-6","IGKV1-9","IGKV1D-39","IGKV2-30","IGKV2D-29","IGKV3-11","IGKV3-15","IGKV3-20","IGKV4-1","IGLC1","IGLC2","IGLC3","IGLJ2","IGLJ3","IGLL5","IGLV1-44","IGLV1-47","IGLV1-51","IGLV2-11","IGLV2-14","IGLV2-23","IGLV2-8","IGLV3-1","IGLV3-19","IGLV3-21","IGLV3-25","IGLV3-9","IGLV4-69","IGLV6-57","JCHAIN","KRT17","LMOD1","LYPD3","MYH11","MYL9","PDLIM3","POSTN","PTGDS","S100P","SFRP2","SFRP4","SNCG","THBS2","TMEM40","TSPAN1"]

    lumInf_bottomn_branch = ["AC132217.4","ACKR1","ACTC1","ADIRF","AKR1C2","C10orf99","CAPNS2","CCL14","CST6","CTB-102L5.9","CXCL17","CYP24A1","DHRS2","FXYD3","HSPB8","IGF2","IGHJ1","IGHV1-18","IGHV1-24","IGHV1-69","IGHV2-5","IGHV3-20","IGHV3-7","IGHV4-31","IGHV4-4","IGHV4-59","IGHV6-1","IGKJ3","IGKJ4","IGKV1-16","IGKV1-17","IGKV1-27","IGKV2-24","IGKV2-28","IGKV3D-15","IGLV1-40","IGLV3-10","IGLV7-46","IGLV8-61","IVL","KRT13","KRT23","KRT7","LY6D","MAL","MMP11","NCCRP1","PADI3","PCP4","PLA2G2A","PSCA","S100A2","S100A9","SYT8","TMEM265","UPK1A","UPK1B","UPK2","UPK3A","UPK3B"]

    lumInf_significant_scatter_plot = ["DES","IGKV1-33","IGLV3-21","KRT23","LY6D","SFRP2","UPK1B"]

    mixed = ["ASRGL1","AURKB","B4GALNT4","B4GALT6","BUB1B","C12orf75","CBX2","CCNE2","CDC25C","CDKAL1","CDKN2A","CENPA","CENPF","CNIH2","COCH","CRABP1","CTD-3193O13.11","DEPDC1","DLGAP5","DLX5","E2F1","E2F3","E2F5","ENHO","ESPL1","EZH2","FAM111B","FOXM1","FSD1","GAL","GGH","GPC2","GTSE1","HES6","HMGB2","HOXB9","KIF15","KIF18B","KIF1A","LINC00668","LMNB1","MELK","MEX3A","MSI1","MT3","MYBL2","MYEF2","NDC80","NEK2","NUF2","NUSAP1","ODC1","PDXP","PIMREG","PLEKHG4B","PLK1","PRAME","RAC3","RAD51AP1","RCOR2","RP11-613M10.6","SBK1","SEPTIN3","SLC29A4","SOX2","STMN1","TMSB15A","TPX2","TUBB2B","UGT8","YBX2","YPEL1"]

    custom_traces.append({"genes": lumInf_top_branch, "title": "LumInf/NS top branch"})
    custom_traces.append({"genes": lumInf_middle_branch, "title": "LumInf/NS middle branch"})
    custom_traces.append({"genes": lumInf_bottomn_branch, "title": "LumInf/NS bottom branch"})
    custom_traces.append({"genes": lumInf_significant_scatter_plot, "title": "LumInf/NS significant scatter"})

    custom_traces.append({"genes": mixed, "title": "Mixed"})

    return custom_traces

def significant_genes(custom_traces):

    small = ['IGKJ1', 'IGKV3-20', 'IGLC2', 'IGKV2-28', 'IGKJ2', 'IGHJ4', 'IGLC1', 'IGKV4-1', 'H19', 'CD14', 'IGHG1', 'C1QB', 'IGHJ5', 'POSTN', 'HTRA3', 'C1QC', 'IGKV1D-39', 'CTHRC1', 'IGHG2', 'AEBP1', 'IGHG3', 'HLA-DQA1', 'IGLC3', 'NNMT', 'CCN2', 'IGKJ4', 'MMP11', 'CTSK']

    # large = ['KRT5', 'DSP', 'S100A2', 'H19', 'LY6D', 'CDH3', 'AQP3', 'MALL', 'KRT16', 'S100A8', 'SERPINB5', 'ANXA8L1', 'SLPI', 'COL17A1', 'KRT14', 'TNS4', 'S100A9']

    large_2 = ['AKR1C2', 'ABO', 'AQP3', 'ANXA8L1', 'ALDH3B2', 'AKR1C3', 'AKR1C1',
       'ADGRF1', 'AATBC', 'AGR2', 'ADGRF4', 'AC124789.1', 'AC003099.2',
       'ADGRV1', 'ADORA2B', 'ANK3', 'ARL4D', 'AREG', 'AC007283.5',
       'ABLIM3', 'ALDH1L1', 'CA12', 'COL17A1', 'GPX2', 'BICDL2', 'ABCA1',
       'CASP1', 'KRT14', 'KRT13', 'ACER2', 'AC011288.2', 'ATP2C2',
       'ALOX5AP', 'ANO1', 'CSTA', 'KRT5', 'AMDHD1', 'ANXA3', 'CA9',
       'CDH3', 'FAM83A', 'CKMT1A', 'C1orf116', 'CELSR2', 'GJB6', 'CKMT1B',
       'ADCY7', 'DEFB1', 'CLDN1', 'AK4', 'LY6D', 'ABI3', 'BATF2',
       'ADAMTSL4', 'AC000123.4', 'AIF1', 'ANKRD18B', 'KRT16', 'FGFR3',
       'CH17-360D5.2', 'CA2', 'AC005537.2', 'C15orf48', 'AC007326.11',
       'ACOT11', 'ALDH1A3', 'ARHGEF35', 'CCL20', 'GJB2', 'ANKRD50',
       'ADSS1', 'ABCC3', 'ARTN', 'AMN', 'ALDH5A1', 'CXCL17', 'ABAT',
       'BNIPL', 'KRT15', 'DSP', 'ANKRD22', 'AC074286.1', 'CLEC2B',
       'CXCL1', 'AFAP1L2', 'CFB', 'ADA', 'GRHL3', 'CAV1', 'AC022154.7',
       'AC004237.1', 'CXCL8', 'CNFN', 'ARNTL2', 'CXCL10', 'CELSR1',
       'IL20RB', 'ADGRG6', 'ALDH1A2', 'LCN2', 'ADAM28', 'GJB5', 'ANGPTL4',
       'ARHGEF37', 'GJB3', 'CFH', 'IVL', 'FAM110C', 'DSC2', 'CCL5',
       'C1orf210', 'C1QB', 'ARHGAP9', 'AC015849.19', 'AP1S3', 'C1orf74',
       'CARD6']
    lumP = ['S100A2', 'PSCA', 'H19', 'LY6D', 'UPK2', 'AQP3', 'ABCC3', 'KRT13', 'IGLC2', 'GPX2', 'CDC42EP5', 'VSIG2', 'DHRS2', 'SPINK1', 'TMPRSS4']

    mixed = ['IGHG1', 'IGLC2', 'CDKN2A']

    lumInf = ['IGKV3-20', 'IGLC2', 'ACTG2', 'IGHJ4', 'IGLC1', 'S100A9', 'H19', 'MUC20', 'IGHG1', 'C1QB', 'POSTN', 'GPX2', 'IGHA1', 'PSCA', 'UPK2', 'AEBP1', 'IGLC3', 'NNMT', 'KRT23', 'CTSK']

    custom_traces.append({"genes":large_2, "title": "Highest large Ba/Sq"})
    custom_traces.append({"genes":small, "title": "Highest small Ba/Sq"})

    custom_traces.append({"genes":lumP, "title": "Highest lump"})
    custom_traces.append({"genes":mixed, "title": "Highest mixed"})
    custom_traces.append({"genes":lumInf, "title": "Highest lumInf"})

    return custom_traces

def low_significant_genes(custom_traces):

    small = ['GSTM2', 'FABP4', 'HSD17B1', 'SCIN', 'NEB', 'AC132217.4', 'S100A3', 'RP11-20D14.6', 'SYNM', 'MYLPF', 'SOX2-OT', 'HES6', 'ALDH1A1', 'CRMP1', 'STAG3', 'CES1', 'LGALS4', 'CCL20', 'FCRLB', 'MSX1', 'SCUBE2', 'TSPAN7', 'IGF2', 'RP11-87N24.3', 'CD79B', 'CTC-425F1.4', 'ITM2A', 'AIFM3', 'NRN1', 'RP5-940J5.9', 'TGM1', 'AKR1C1', 'FABP3', 'RASD1', 'IL20RB', 'SPINK1', 'MMP9', 'SELL', 'GLDC', 'TESC', 'LINC01088', 'CD36']

    large = ['AC005301.9', 'S100A3', 'RP11-20D14.6', 'SYNM', 'MYLPF', 'SOX2-OT', 'HES6', 'ALDH1A1', 'STAG3', 'HSPB6', 'CES1', 'LGALS4', 'CCL20', 'FCRLB', 'MSX1', 'TSPAN7', 'EFHD1', 'RP11-87N24.3', 'CD79B', 'CTC-425F1.4', 'ITM2A', 'AIFM3', 'RP5-940J5.9', 'FABP3', 'RASD1', 'MMP9', 'APLP1', 'SELL', 'ID4', 'SERPINI1', 'TESC', 'LINC01088', 'CD36']

    large_2 = ['AEBP1', 'AK4', 'ADAM19', 'ALOX5AP', 'ADA', 'ANXA3', 'ADCY7',
       'AIF1', 'AREG', 'ADAMTSL4', 'ADAMTS2', 'KRT14', 'C1QB', 'ACTG2',
       'C1QC', 'CAV1', 'KRT5', 'ARNTL2', 'A1BG', 'CDH3', 'ANKH',
       'ADORA2B', 'ARL4D', 'BCL2A1', 'CCDC80', 'ANXA8L1', 'AFAP1L2',
       'CD14', 'CXCL10', 'CHI3L1', 'ALDH1A1', 'CCL4', 'CDA', 'KRT16',
       'ABI3', 'BIN1', 'CHST11', 'CXCL1', 'CASP1', 'CCL20', 'AGPAT4',
       'CTHRC1', 'AMOTL1', 'CHST15', 'CCL3', 'CA9', 'CD163', 'DSC2', 'F3',
       'C15orf48', 'COL17A1', 'BATF2', 'DSP', 'ARHGAP9', 'EFEMP1', 'AXL',
       'COL5A1', 'CCL5', 'ANPEP', 'ALDH2', 'CCDC71L', 'CXCL8', 'C12orf75',
       'COL5A2', 'CTSV', 'DEFB1', 'CD53', 'GJB2', 'COL16A1', 'IL20RB',
       'CES1', 'AC018816.3', 'FAM83A', 'COPZ2', 'CLEC2B', 'CA2',
       'AGAP2-AS1', 'C1orf162', 'CFB', 'ALDH1A3', 'ADAP2', 'AOAH',
       'ADCY9', 'AMDHD1', 'BDKRB2', 'ABCA1', 'CALHM6', 'FCGR3A',
       'CFAP251', 'ANO1', 'CD7', 'APOC2', 'AC007283.5', 'ANGPTL2', 'DPYD',
       'FAP', 'AQP3', 'ARHGEF37', 'IDO1', 'ADAMTS4', 'CD86', 'BIRC3',
       'CRYAB', 'CD300A', 'CFI', 'ARHGAP22', 'CELSR2', 'CELSR1',
       'AFAP1L1', 'CSF1R', 'CDCA7L', 'ATP10D', 'COL5A3', 'FCGR2A', 'ASPM',
       'CD3E', 'C1QTNF1', 'CCN2', 'C1orf74', 'ECM1', 'CD2', 'LAMA3',
       'CD52', 'GNLY', 'HTRA3', 'CYBB', 'FAM20A', 'CLCF1', 'CSF1',
       'LAMC2', 'CEACAM19', 'CCL2', 'CGAS', 'IGHG1', 'DSE', 'EPSTI1']


    lumP =  ['MT1M', 'HSD17B1', 'CD36', 'PTHLH', 'LMOD1', 'KRT14', 'S100A3', 'RP11-20D14.6', 'SYNM', 'CHI3L1', 'SOX2-OT', 'ALDH1A1', 'SUGCT', 'IL1R2', 'IGHG4', 'STAG3', 'OLR1', 'CES1', 'SPP1', 'CCL20', 'TSPAN7', 'EFHD1', 'GZMB', 'ANPEP', 'UCHL1', 'COMP', 'CD79B', 'ITM2A', 'AREG', 'TGM1', 'FABP3', 'EFEMP1', 'RASD1', 'IL20RB', 'MMP9', 'APLP1', 'SELL', 'IDO1', 'IGHM']

    mixed = ['MT1M', 'MMP1', 'COL17A1', 'PTHLH', 'LMOD1', 'KRT14', 'S100A3', 'SYNM', 'SOX2-OT', 'SUGCT', 'IL1R2', 'MMP7', 'CRMP1', 'OLR1', 'HSPB6', 'CES1', 'CCL20', 'CNN1', 'RP11-54H7.4', 'MYH11', 'GABRP', 'RP11-87N24.3', 'RRAD', 'CACNA2D4', 'CD79B', 'CTC-425F1.4', 'ITM2A', 'RP5-940J5.9', 'TGM1', 'FABP3', 'IL20RB', 'MMP9', 'LINC01088', 'CD36']

    lumInf = ['AC005301.9', 'DPY19L2', 'LETM2', 'CD36', 'PTHLH', 'KRT14', 'S100A3', 'RP11-20D14.6', 'RP3-467D16.3', 'SOX2-OT', 'HES6', 'ALDH1A1', 'CRMP1', 'STAG3', 'CES1', 'CCL20', 'SNORC', 'MSX1', 'RP11-54H7.4', 'RP11-87N24.3', 'CDK5R1', 'RRAD', 'CD79B', 'CTC-425F1.4', 'ITM2A', 'GLDC', 'RP5-940J5.9', 'TGM1', 'FABP3', 'RASD1', 'IL20RB', 'APLP1', 'ID4', 'IDO1', 'LINC01088', 'SCGB3A2']

    # There are a lot of genes that have low expression and are shared across. So, we're interested in the genes that are particulary low for a subgroup
    all_genes = set(small) | set(large) | set(lumP) | set(mixed) | set(lumInf) 

    custom_traces.append({"genes": set(large_2), "title": "Lowest arge Ba/Sq"})
    custom_traces.append({"genes": all_genes & set(small), "title": "Lowest small Ba/Sq"})

    custom_traces.append({"genes": all_genes & set(lumP), "title": "Lowest lump"})
    custom_traces.append({"genes": all_genes & set(mixed), "title": "Lowest mixed"})
    custom_traces.append({"genes": all_genes & set(lumInf), "title": "Lowest lumInf"})

    return custom_traces


from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from scanpy.pl._utils import savefig_or_show

"""
Dev notes: Preliminary plotting scripts -- needs testing, show/save implementation
split_violins now runs on the tutorial data and on my data 
"""
class DialoguePlot:
    """Plotting functions for Dialogue."""

    @staticmethod
    def split_violins(
        adata: AnnData, 
        celltype_key: str, 
        split: str, 
        mcp: str="mcp_0", 
        split_tuple: (str,str) | None = None,
        title: str |None = None,
        show: bool | None = None,
        save: bool | str | None = None,
        return_fig: bool = False,
        **kwargs,
        ):
        """Plot MCP results from Dialogue on split violins, with one violin per cell type.

        Args:
            adata: AnnData with MCP scores
            celltype_key: the key in AnnData.obs used to create the violins
            split: the key in AnnData.obs used to split the violins. Must have exactly 2 values if split_tuple is None
            mcp: the key in AnnData.obs which contains the MCP score being plotted, in the form "mcp_x"
            groupby: the key in AnnData.obs which contains the groups of interest (typically cell type)
            split_tuple: None (if AnnData.obs['split'] has exactly 2 values) or a tuple of strings to subset to only
                cells with those strings in AnnData.obs['split']
            title:
            show:Show the plot, do not return axis.
            save: If `True` or a `str`, save the figure. A string is appended to the default filename.
                  Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
            **kwargs: Additional arguments passed to `sns.violinplot`

        """
        if split_tuple is not None:
            df = sc.get.obs_df(adata[adata.obs[split].isin(split_tuple)], [celltype_key,mcp,split])
        else:
            df = sc.get.obs_df(adata, [celltype_key,mcp,split])

        # work save/show into it 
        p2 = sns.violinplot(data=df, x=celltype_key,y=mcp,hue=split, split=True, **kwargs)
        p2.set_xticklabels(p2.get_xticklabels(), rotation=90)
        p2.set(title=title)
        savefig_or_show("split_violin_"+mcp+"_", show=show, save=save)
        if return_fig:
            return(p2)


    @staticmethod
    def mcp_pairplot(
        adata: AnnData,
        celltype_key: str ,
        sample_key: str,
        mcp: str="mcp_0",
        title: str | None = None,
        color: str | None = None,
        show: bool | None = None,
        save: bool | str | None = None,
        **kwargs,
        ):
        """Plot MCP results from Dialogue on a pair plot as in DIALOGUE 

        Args:
            adata: AnnData with MCP scores
            celltype_key: the key in AnnData.obs which contains the celltype key 
            color: the key in AnnData.obs used to color the points and histograms (optional)
            mcp: the key in AnnData.obs which contains the MCP score being plotted
            groupby: the key in AnnData.obs which contains the groups of interest (typically cell type)
            split: the key in AnnData.obs which contains the info being used to split 
            title:
            show:Show the plot, do not return axis.
            save: If `True` or a `str`, save the figure. A string is appended to the default filename.
                  Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
            **kwargs: Additional arguments passed to `sns.pairplot`

        """

        ## TODO: figure out how to move the legend to the main plot
        ## OR: print the Pearson correlations like DLG does


        mean_mcp = adata.obs.groupby([sample_key,celltype_key])[mcp].mean()
        mean_mcp = mean_mcp.reset_index()
        mcp_pivot = pd.pivot(mean_mcp[[sample_key,celltype_key,mcp]],index=sample_key,columns = celltype_key )[mcp]
        if color is not None:
            mcp_pivot = pd.concat([mcp_pivot,adata.obs.groupby([sample_key]).agg(pd.Series.mode)[color]],axis=1)
            p2 = sns.pairplot(mcp_pivot, corner=True, hue=color, **kwargs)

        else: 
            p2 = sns.pairplot(mcp_pivot, corner=True, **kwargs)


        if title is not None:
            p2.fig.suptitle(title)
        else:
            p2.fig.suptitle("MCP score "+mcp)

        savefig_or_show("MCP_pairplot_"+mcp+"_", show=show, save=save)




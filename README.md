# Endo-atlas

This repo contains the code used to construct the endo-atlas dash app. Assets and obj folders were ignored due to file size constraint.

The app is available at http://endoatlastest-env.eba-qtpfn7qz.us-east-1.elasticbeanstalk.com/ .


## Overview

The app was built with `dash`.  A simplified 3D model of Nematostella embryo was constructed with `blender®`. The expression pattern of 51 marker genes was summarized to binarized values at 8993 locations along the surface of this 3D model. The spatial expression patterns of other genes were then reconstructed with `NovoSpaRc` (Nitzan et al., 2019) using single cell RNAseq data. The predicted expression patterns were then visualized in 3D using `plotly`.

The app allows the user to select genes of interest, the outputs include in silico prediction of expression in 3D, available in situ image, and top 10 genes with similar spatial pattern based on pattern correlation analysis.


<p align="center">
  <img src="https://github.com/wanqingshao/endo-atlas/blob/main/app_overview.png?raw=true" width="1000" alt="app overview">
</p>

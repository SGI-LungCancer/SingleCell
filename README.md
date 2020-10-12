# Single-cell RNA sequencing for lung adenocarcinoma

In order to separate maligant tumor cells from non-malignant cells, we calculate CNV aberrations inferring the perturbation of chromosomal gene expression.

1. Adjusting the proportion of putative malignant cells below 20%.
2. Filter out less infermative genes (default : less than 10 cells and mean expression of less than 0.1 at log2 scale).
3. Transformation into Z-score and limit the scale -3 to 3.
4. Sorting the genes by their chromosomal position and estimage CNV signals using the window size (default = 100 genes).
5. Summerize CNV signal with two parameters and classify malignant cells and non-malignant cells.
  - CNV signals (MS, Mean of Squares) : Mean squares of estimates across all windows.
  - CORR (Correlation with the high CNV signal cells) : Correlation of the CNV of each cell with the average of the top5% cells.
  - Maligant cells were classified if their CNV signals (MS) > 0.02 or CNV correlation (CORR) > 0.2.
  
  
# Example code

After download the codes, run example data using 'calc_Chromosomal_Expression_Pattern.R'.




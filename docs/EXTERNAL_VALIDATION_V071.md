# External validation v0.7.1

## Completed S-BIAD3159 run

The locked four-field workflow completed on `Upper_Left`, `Upper_Right`, `Bottom_Left`, and
`Bottom_Right`, with red-channel correction α ∈ {0, 0.5, 1}. All four fields preserved a negative
association between relative Hoechst intensity and distance to structural CD31 signal, and all four
preserved enrichment within 10 µm relative to 50–100 µm. Across the 12 sensitivity runs, median
Spearman was **-0.2049** and median near/far
ratio was **3.062×**.

This validates the relative extraction method in independent LLC mouse-tumour fields. It does not
identify perfused vessels in RCC and does not measure a therapeutic antibody.

## Native breast Xenium ERBB2 extraction

The provider `cell_features` Zarr schema was decoded using group attributes (`feature_keys`,
`number_features`, `number_cells`), packed integer cell IDs, and feature-by-cell CSR with CSC fallback.
The two sections contributed **679,197** labelled cells. Provider-labelled tumour cells had mean ERBB2
RNA 2.7289 in HER2-2+ and 56.2783 in HER2-3+, a descriptive **20.623×** ratio. Tumour-cell positive
fractions were 0.7769 and 0.9878, respectively.

The comparison is between two independent sections. Cells are not biological replicates; no
cell-level inferential p-value is reported. ERBB2 RNA is not a surface-receptor count and the source Cy5
calibration is not transferred.

## Aggregate status

`EXTERNAL_PERFUSION_ALL_FOUR_FIELDS_AND_BREAST_XENIUM_ERBB2_RNA_VALIDATED_MODEL_PHARMACOLOGICAL_CONCORDANCE_NOT_COMPUTED`

Absolute RCC outputs remain `NOT_COMPUTED`.

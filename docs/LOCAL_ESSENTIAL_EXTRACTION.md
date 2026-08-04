# Low-memory extraction on a local computer

The 36.1 GB Xenium ZIP does **not** need to be uploaded or loaded into RAM. The standalone extractor uses
only the Python standard library, reads the ZIP central directory, and streams only the cell-level files
needed by `reach-gap`.

## Requirements

- Python 3.11 or newer.
- The original Xenium ZIP stored on a local disk.
- Free disk space for the selected outputs. The script reports the exact selected uncompressed size with
  `--list-only` before extracting.
- No GPU and no scientific Python environment are required.

## Windows PowerShell

```powershell
python .\extract_xenium_essentials_standalone.py `
  "D:\Xenium\Xenium_V1_Human_Kidney_FFPE_Protein_updated_outs.zip" `
  --output-dir "D:\Xenium\xenium-essential-package" `
  --list-only

python .\extract_xenium_essentials_standalone.py `
  "D:\Xenium\Xenium_V1_Human_Kidney_FFPE_Protein_updated_outs.zip" `
  --output-dir "D:\Xenium\xenium-essential-package" `
  --part-size-mb 95
```

## macOS or Linux

```bash
python3 extract_xenium_essentials_standalone.py \
  /path/Xenium_V1_Human_Kidney_FFPE_Protein_updated_outs.zip \
  --output-dir /path/xenium-essential-package \
  --list-only

python3 extract_xenium_essentials_standalone.py \
  /path/Xenium_V1_Human_Kidney_FFPE_Protein_updated_outs.zip \
  --output-dir /path/xenium-essential-package \
  --part-size-mb 95
```

## Memory model

Only an 8 MiB data chunk plus ZIP decompressor buffers are held at once. A large HDF5 member is written
directly as 95 MB parts; no temporary full HDF5 copy is created. The process is therefore expected to
fit comfortably within 16 GB RAM. Disk throughput, source ZIP integrity and free output storage are the
practical constraints.

## Upload after completion

Upload the complete `xenium-essential-package/` directory to Drive. It should include:

- `essential_package_manifest.json`;
- `selected_essential_members.json`;
- `xenium_bundle_inventory.json`;
- `files/` containing small selected members;
- `split/cell_feature_matrix.h5/` and its split manifest when the HDF5 exceeds 95 MB;
- any split cell table and its manifest if applicable.

Do not rename parts or omit the manifests. The full 36.1 GB source ZIP is not needed after this gate.

quick filtering pipeline written with extensive agent usage.

The remainder of the readme is agent-written, not necessarily endorsed.


# filter_h5ads

A robust, type-safe pipeline for filtering h5ad files from CRISPR perturbation screens.

## Features

- **Immutable Configuration**: Pydantic-based configs ensure reproducibility
- **Deterministic Naming**: Config hashing creates unique, traceable output files
- **Comprehensive Filtering**: UMI counts, guide quality, mitochondrial content, gene detection
- **Metadata Filtering**: Include/exclude cells by `.obs` values (e.g., cell type)
- **Rich Inspection**: Detailed overview of h5ad file contents
- **Batch Processing**: Process multiple files with the same pipeline
- **Full Provenance**: Filtering metadata stored in output files
- **Type Safety**: Strict typing for Python 3.13+

## Installation

```bash
pip install -r requirements.txt
```

This will install the package and make the `filter-h5ads` CLI command available.

## Quick Start

### Command Line Interface

Inspect an h5ad file:
```bash
filter-h5ads inspect data.h5ad
filter-h5ads inspect data.h5ad --n-examples 10
filter-h5ads inspect data.h5ad --backed  # Memory efficient backed mode
```

Filter a single file:
```bash
filter-h5ads filter data.h5ad config.json
filter-h5ads filter data.h5ad config.json --output-dir filtered/
filter-h5ads filter data.h5ad config.json --backed  # Memory efficient
```

Batch process multiple files:
```bash
filter-h5ads batch data/raw/ config.json --output-dir data/filtered/
filter-h5ads batch data/raw/ config.json --output-dir data/filtered/ --backed
```

### Backed Mode (Memory Efficient)

For large h5ad files, you can use backed mode to reduce memory usage:

**CLI**: Add `--backed` flag to any command
```bash
filter-h5ads inspect huge_file.h5ad --backed
```

**Python API**: Pass `backed=True` to `load_h5ad()`
```python
adata = load_h5ad("huge_file.h5ad", backed=True)
```

**Environment Variable**: Set globally for all operations
```bash
export LOAD_AD_BACKED=1
filter-h5ads inspect data.h5ad  # Will use backed mode
```

### Python API

```python
from pathlib import Path
from filter_h5ads import (
    load_h5ad,
    print_h5ad_overview,
    FilterPipelineConfig,
    GeneDetectionFilterConfig,
    GroupCountFilterConfig,
    GroupCountAnnotateConfig,
    ObsColumnTransformConfig,
    PseudobulkConfig,
    UMIFilterConfig,
    GuideFilterConfig,
    MitochondrialFilterConfig,
    ObsValueFilterConfig,
    run_pipeline,
    save_filtered_h5ad,
)

# Load and inspect
adata = load_h5ad("data.h5ad")
print_h5ad_overview(adata)

# Configure pipeline
config = FilterPipelineConfig(
    pipeline_name="standard_crispr_qc",
    obs_column_transform=ObsColumnTransformConfig(
        input_columns=["drugname_drugconc"],
        operations=[
            {"op": "parse_python_literal"},
            {"op": "first"},
            {"op": "pick_indices", "indices": [1, 2]},
            {"op": "join", "sep": " "},
        ],
        output_column="dose",
    ),
    obs_value_filter=ObsValueFilterConfig(key="cell-type", values=["doublet"], exclude=True),
    umi_filter=UMIFilterConfig(min_counts=15000),
    guide_filter=GuideFilterConfig(guide_column='pass_guide_filter'),
    mito_filter=MitochondrialFilterConfig(max_pct_mt=20.0),
    gene_filter=GeneDetectionFilterConfig(min_genes=200),
    group_count_annotate=GroupCountAnnotateConfig(
        groupby=["cell_line", "drug", "dose"],
        output_column="n_cells_in_condition",
    ),
    group_count_filter=GroupCountFilterConfig(
        groupby=["cell_line", "drug", "dose"],
        min_cells=50,
    ),
    pseudobulk=PseudobulkConfig(
        groupby=["cell_line", "drug", "dose"],
        obs_aggregations={
            "batch": {"agg": "assert_constant"},
            "total_counts": {"agg": "sum", "output_column": "total_counts_sum"},
            "pct_counts_mt": {"agg": "mean"},
            "cell_barcode": {"agg": "list"},
            "gene_target": {"agg": "proportions", "output_column": "gene_target_props"},
        },
    ),
)

# Run pipeline
adata_filtered, stats = run_pipeline(adata, config)

# Save with automatic naming: data_filtered_{hash}.h5ad
output_path = save_filtered_h5ad(adata_filtered, Path("data.h5ad"), config)
```

## Pipeline Configuration

Each filter can be independently configured and enabled/disabled. The pipeline executes in order: **Ensembl Conversion** → **Obs Column Transform** → **Obs Value** → UMI → Guide → Mito → Gene → **Group Count Annotate** → **Group Count Filter** → **Pseudobulk**.

### Ensembl ID Conversion (Preprocessing)
```python
from comlm.config.data_configs import SafeGeneSymbolRemappingRef

EnsemblConversionConfig(
    gene_remapping_ref=SafeGeneSymbolRemappingRef.get_default_for_vcc(),
    mask_getter=lambda: get_allowed_genes(),  # Your function to get allowed genes
    source_target_column="gene_target",        # Column to convert
    target_column_name="target_gene",          # New column name
    enabled=True
)
```

This step:
1. Creates a new `target_gene` column with Ensembl IDs (converts "Non-Targeting" → "non-targeting")
2. Filters cells to only keep those with targets in the allowed gene set
3. Converts `.var_names` to Ensembl IDs
4. Adds zero columns for missing genes in the mask
5. Ensures all genes in mask are present in the data

### Obs Value Filter (e.g., Cell Type)
```python
ObsValueFilterConfig(
    key="cell-type",                 # Column name in .obs
    values=["doublet", "low_qc"],     # Values to match against
    exclude=True,                    # True=drop matches, False=keep only matches
    enabled=True,
)
```

### Obs Column Transform (derive new metadata)
```python
ObsColumnTransformConfig(
    input_columns=["drugname_drugconc"],
    operations=[
        {"op": "parse_python_literal"},            # if values are stored as strings
        {"op": "first"},                           # [('Drug', 5.0, 'uM')] -> ('Drug', 5.0, 'uM')
        {"op": "pick_indices", "indices": [1, 2]}, # -> (5.0, 'uM')
        {"op": "join", "sep": " "},                # -> "5.0 uM"
    ],
    output_column="dose",
    overwrite=False,
    on_error="raise",
)
```

### Pseudobulk (group + sum counts)
```python
PseudobulkConfig(
    groupby=["cell_line", "drug", "dose"],  # any unique combination becomes 1 pooled sample
    layer=None,                              # sum `.X` (set to a layer name to sum `.layers[layer]`)
    obs_aggregations={
        # For each non-groupby `.obs` column you care about, specify how to aggregate it:
        "batch": {"agg": "assert_constant"},
        "total_counts": {"agg": "sum", "output_column": "total_counts_sum"},
        "pct_counts_mt": {"agg": "mean"},
        "cell_barcode": {"agg": "list"},  # list of values across pooled cells
        "gene_target": {"agg": "proportions", "output_column": "gene_target_props"},
        # Unspecified `.obs` columns are dropped by default.
    },
)
```

### Group Count Filter (min cells per condition)
```python
GroupCountFilterConfig(
    groupby=["cell_line", "drug", "dosage"],
    min_cells=50,  # keep only conditions with ≥50 cells
)
```

### Group Count Annotate (group size per cell)
```python
GroupCountAnnotateConfig(
    groupby=["cell_line", "drug", "dosage"],
    output_column="n_cells_in_condition",
)
```

### UMI Count Filter
```python
UMIFilterConfig(
    min_counts=15000,           # Minimum UMI threshold
    count_column="total_counts", # Column name
    enabled=True                 # Can disable
)
```

### Guide Quality Filter
```python
GuideFilterConfig(
    guide_column="pass_guide_filter",  # Column name
    required_value=1,                   # Required value (1 or True)
    enabled=True
)
```

### Mitochondrial Content Filter
```python
MitochondrialFilterConfig(
    max_pct_mt=20.0,              # Max % mitochondrial
    mt_column="pct_counts_mt",    # Column name
    enabled=True
)
```

### Gene Detection Filter
```python
GeneDetectionFilterConfig(
    min_genes=200,                    # Min genes detected
    gene_column="n_genes_by_counts",  # Column name
    enabled=True
)
```

## Batch Processing

Process multiple files with the same configuration:

```python
from filter_h5ads import batch_filter_h5ads

results = batch_filter_h5ads(
    input_paths=Path("data/raw"),  # Directory or list of files
    config=config,
    output_dir=Path("data/filtered"),
    skip_existing=True  # Skip already-processed files
)
```

## Output Naming Convention

Files are named based on config hash for reproducibility:

```
original.h5ad → original_filtered_a3f9c2e1.h5ad
              → original_filtered_a3f9c2e1.json (config)
```

Different configurations produce different hashes:
- Same config always produces same hash (reproducible)
- Different configs produce different hashes (traceable)

## Inspection Tools

### Print Overview
```python
print_h5ad_overview(adata, n_examples=5)
```

Shows:
- Dimensions (cells × genes)
- All `.obs` columns with types and examples
- All `.var` columns with types and examples
- Layers, embeddings, and other data structures

### Get Programmatic Summary
```python
summary = get_h5ad_summary(adata)
# Returns structured dict with all information
```

### Print Filtering Results
```python
print_filtering_summary(stats)
```

Shows:
- Cells retained at each step
- Breakdown by `gene_target`
- Overall pipeline statistics

## Provenance Tracking

All filtering operations are recorded in `adata.uns['filtering_provenance']`:

```python
{
    "pipeline_name": "standard_crispr_qc",
    "pipeline_hash": "a3f9c2e1",
    "timestamp": "2025-10-16T10:30:00",
    "config": {...},
    "filter_stats": [...],
    "initial_cells": 250000,
    "final_cells": 198450
}
```

## Architecture

```
filter_h5ads/
├── config.py       # Pydantic configuration models
├── filters.py      # Individual filter implementations
├── pipeline.py     # Pipeline orchestration
├── io.py          # File loading/saving/batch processing
├── inspection.py  # Overview and summary tools
└── utils.py       # Helper functions
```

## Examples

See [example_usage.py](example_usage.py) for comprehensive examples:

1. Single file filtering
2. Batch processing
3. Custom pipelines
4. Selectively disabling filters
5. Multiple configurations

## Type Safety

All code uses strict typing compatible with Python 3.13:
- Native type hints (`list`, `dict`, `tuple`)
- Pydantic for validation
- TYPE_CHECKING for circular imports
- Full mypy/pyright compatibility

## Dependencies

- `anndata` >= 0.10.0 - Core data structure
- `pydantic` >= 2.0.0 - Configuration management
- `numpy` >= 1.24.0 - Numerical operations
- `pandas` >= 2.0.0 - DataFrame operations

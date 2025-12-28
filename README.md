quick filtering pipeline written with extensive agent usage.

The remainder of the readme is agent-written, not necessarily endorsed.


# filter_h5ads

A robust, type-safe pipeline for filtering h5ad files from CRISPR perturbation screens.

## Features

- **Immutable Configuration**: Pydantic-based configs ensure reproducibility
- **Deterministic Naming**: Config hashing creates unique, traceable output files
- **Comprehensive Filtering**: UMI counts, guide quality, mitochondrial content, gene detection
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
    UMIFilterConfig,
    GuideFilterConfig,
    run_pipeline,
    save_filtered_h5ad,
)

# Load and inspect
adata = load_h5ad("data.h5ad")
print_h5ad_overview(adata)

# Configure pipeline
config = FilterPipelineConfig(
    pipeline_name="standard_crispr_qc",
    umi_filter=UMIFilterConfig(min_counts=15000),
    guide_filter=GuideFilterConfig(guide_column='pass_guide_filter'),
    mito_filter=MitochondrialFilterConfig(max_pct_mt=20.0),
    gene_filter=GeneDetectionFilterConfig(min_genes=200)
)

# Run pipeline
adata_filtered, stats = run_pipeline(adata, config)

# Save with automatic naming: data_filtered_{hash}.h5ad
output_path = save_filtered_h5ad(adata_filtered, Path("data.h5ad"), config)
```

## Pipeline Configuration

Each filter can be independently configured and enabled/disabled. The pipeline executes in order: **Ensembl Conversion** → UMI → Guide → Mito → Gene.

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

"""I/O operations for loading, saving, and batch processing h5ad files.

Handles file operations with proper naming conventions based on config hashing.
"""

import os
from pathlib import Path
from typing import Any

from anndata import AnnData, read_h5ad
from loguru import logger

from filter_h5ads.config import FilterPipelineConfig
from filter_h5ads.pipeline import run_pipeline, validate_pipeline_config

# Global flag to control backed mode loading
# Set LOAD_AD_BACKED=1 environment variable to load h5ad files in backed mode
LOAD_AD_BACKED = os.environ.get("LOAD_AD_BACKED", "0") == "1"


def load_h5ad(filepath: str | Path, backed: bool | None = None) -> AnnData:
    """Load an h5ad file with validation.

    Args:
        filepath: Path to the h5ad file
        backed: If True, load in backed mode. If None, uses LOAD_AD_BACKED env var.

    Returns:
        Loaded AnnData object

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file is not a valid h5ad

    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    if not filepath.suffix == ".h5ad":
        raise ValueError(f"File must have .h5ad extension: {filepath}")

    # Determine backed mode: explicit param > environment variable > default (False)
    use_backed = backed if backed is not None else LOAD_AD_BACKED

    logger.info(f"Loading {filepath}" + (" (backed mode)" if use_backed else ""))

    try:
        if use_backed:
            adata = read_h5ad(filepath, backed="r")
        else:
            adata = read_h5ad(filepath)
        logger.info(f"Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        return adata
    except Exception as e:
        raise ValueError(f"Failed to load h5ad file: {e}") from e


def save_filtered_h5ad(
    adata: AnnData,
    original_path: Path,
    config: FilterPipelineConfig,
    output_dir: Path | None = None,
    compression: str = "gzip",
) -> Path:
    """Save filtered h5ad file with deterministic naming based on config hash.

    Creates files:
    - {name}_filtered_{hash}.h5ad: The filtered data
    - {name}_filtered_{hash}.json: The configuration used

    Args:
        adata: Filtered AnnData object to save
        original_path: Path to the original h5ad file (for naming)
        config: Pipeline configuration that was used
        output_dir: Directory to save outputs (defaults to same as original)
        compression: Compression method for h5ad file

    Returns:
        Path to the saved h5ad file

    """
    original_path = Path(original_path)
    config_hash = config.get_hash()

    # Determine output directory
    if output_dir is None:
        output_dir = original_path.parent
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # Create output filenames
    base_name = original_path.stem
    h5ad_filename = f"{base_name}_filtered_{config_hash}.h5ad"
    config_filename = f"{base_name}_filtered_{config_hash}.json"

    h5ad_path = output_dir / h5ad_filename
    config_path = output_dir / config_filename

    # Save h5ad file
    logger.info(f"Saving filtered h5ad to: {h5ad_path}")
    adata.write_h5ad(h5ad_path, compression=compression)

    # Save configuration
    logger.info(f"Saving configuration to: {config_path}")
    config.save_json(config_path)

    logger.info(f"Successfully saved {adata.n_obs:,} cells")

    return h5ad_path


def batch_filter_h5ads(
    input_paths: list[Path] | Path,
    config: FilterPipelineConfig,
    output_dir: Path,
    skip_existing: bool = True,
    compression: str = "gzip",
    backed: bool | None = None,
) -> list[dict[str, Any]]:
    """Process multiple h5ad files with the same filtering pipeline.

    Args:
        input_paths: List of paths to h5ad files, or a directory containing
            h5ad files
        config: Pipeline configuration to apply to all files
        output_dir: Directory to save filtered outputs
        skip_existing: Skip files that have already been processed (based on
            hash)
        compression: Compression method for h5ad files
        backed: If True, load files in backed mode. If None, uses
            LOAD_AD_BACKED env var.

    Returns:
        List of dictionaries with results for each file processed

    """
    # Handle directory input
    if isinstance(input_paths, Path) and input_paths.is_dir():
        input_paths = sorted(input_paths.glob("*.h5ad"))
    elif isinstance(input_paths, (str, Path)):
        input_paths = [Path(input_paths)]
    else:
        input_paths = [Path(p) for p in input_paths]

    if not input_paths:
        logger.warning("No h5ad files found to process")
        return []

    logger.info(f"Processing {len(input_paths)} h5ad files")
    logger.info(f"Output directory: {output_dir}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results: list[dict[str, Any]] = []
    config_hash = config.get_hash()

    for i, input_path in enumerate(input_paths, 1):
        logger.info(f"\n{'=' * 80}")
        logger.info(f"Processing file {i}/{len(input_paths)}: {input_path.name}")
        logger.info(f"{'=' * 80}")

        # Check if already processed
        expected_output = output_dir / f"{input_path.stem}_filtered_{config_hash}.h5ad"
        if skip_existing and expected_output.exists():
            logger.info(f"Skipping (already exists): {expected_output.name}")
            results.append(
                {
                    "input_path": str(input_path),
                    "output_path": str(expected_output),
                    "status": "skipped",
                    "reason": "already_exists",
                },
            )
            continue

        try:
            # Load data
            adata = load_h5ad(input_path, backed=backed)

            # Validate configuration
            validation_errors = validate_pipeline_config(config, adata)
            if validation_errors:
                error_msg = "; ".join(validation_errors)
                logger.error(f"Validation failed: {error_msg}")
                results.append(
                    {
                        "input_path": str(input_path),
                        "status": "failed",
                        "reason": "validation_error",
                        "errors": validation_errors,
                    },
                )
                continue

            # Run pipeline
            adata_filtered, stats = run_pipeline(adata, config)

            # Save results
            output_path = save_filtered_h5ad(
                adata_filtered,
                input_path,
                config,
                output_dir=output_dir,
                compression=compression,
            )

            # Record success
            results.append(
                {
                    "input_path": str(input_path),
                    "output_path": str(output_path),
                    "status": "success",
                    "initial_cells": stats[0]["n_cells_before"]
                    if stats
                    else adata.n_obs,
                    "final_cells": adata_filtered.n_obs,
                    "filter_stats": stats,
                },
            )

            logger.info(f"✓ Successfully processed: {input_path.name}")

        except Exception as e:
            logger.error(f"✗ Failed to process {input_path.name}: {e}")
            results.append(
                {
                    "input_path": str(input_path),
                    "status": "failed",
                    "reason": "processing_error",
                    "error": str(e),
                },
            )

    # Print summary
    _print_batch_summary(results)

    return results


def _print_batch_summary(results: list[dict[str, Any]]) -> None:
    """Print summary of batch processing results.

    Args:
        results: List of result dictionaries from batch processing

    """
    print("\n" + "=" * 80)
    print("BATCH PROCESSING SUMMARY")
    print("=" * 80)

    total = len(results)
    success = sum(1 for r in results if r["status"] == "success")
    skipped = sum(1 for r in results if r["status"] == "skipped")
    failed = sum(1 for r in results if r["status"] == "failed")

    print(f"\nTotal files: {total}")
    print(f"  ✓ Success: {success}")
    print(f"  ⊘ Skipped: {skipped}")
    print(f"  ✗ Failed:  {failed}")

    if success > 0:
        print("\nSuccessfully processed:")
        for result in results:
            if result["status"] == "success":
                initial = result.get("initial_cells", 0)
                final = result.get("final_cells", 0)
                pct = (final / initial * 100) if initial > 0 else 0
                input_name = Path(result["input_path"]).name
                print(f"  • {input_name}: {initial:,} → {final:,} cells ({pct:.1f}%)")

    if failed > 0:
        print("\nFailed:")
        for result in results:
            if result["status"] == "failed":
                input_name = Path(result["input_path"]).name
                reason = result.get("reason", "unknown")
                print(f"  • {input_name}: {reason}")

    print("=" * 80)


def list_filtered_files(
    directory: Path, original_name: str | None = None
) -> list[Path]:
    """List all filtered h5ad files in a directory.

    Args:
        directory: Directory to search
        original_name: Optional original filename to filter by

    Returns:
        List of paths to filtered h5ad files

    """
    directory = Path(directory)

    if original_name:
        # Look for files matching this original name
        pattern = f"{Path(original_name).stem}_filtered_*.h5ad"
    else:
        # Look for all filtered files
        pattern = "*_filtered_*.h5ad"

    return sorted(directory.glob(pattern))


def get_config_for_filtered_file(h5ad_path: Path) -> FilterPipelineConfig | None:
    """Load the configuration file associated with a filtered h5ad file.

    Args:
        h5ad_path: Path to filtered h5ad file

    Returns:
        FilterPipelineConfig if config file exists, None otherwise

    """
    h5ad_path = Path(h5ad_path)
    config_path = h5ad_path.with_suffix(".json")

    if config_path.exists():
        return FilterPipelineConfig.load_json(config_path)

    return None

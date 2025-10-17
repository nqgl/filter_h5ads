#!/usr/bin/env python3
"""
Example usage of the filter_h5ads package.

Demonstrates how to:
1. Load and inspect an h5ad file
2. Configure a filtering pipeline
3. Run the pipeline on a single file
4. Process multiple files in batch
"""

from pathlib import Path

from filter_h5ads import (
    FilterPipelineConfig,
    GeneDetectionFilterConfig,
    GuideFilterConfig,
    MitochondrialFilterConfig,
    UMIFilterConfig,
    batch_filter_h5ads,
    load_h5ad,
    print_filtering_summary,
    print_h5ad_overview,
    run_pipeline,
    save_filtered_h5ad,
)


def example_single_file():
    """Example: Filter a single h5ad file."""
    print("EXAMPLE 1: Single File Filtering")
    print("=" * 80)

    # Load the h5ad file
    input_path = Path("data/example.h5ad")

    # Check if file exists (for demo purposes)
    if not input_path.exists():
        print(f"Note: {input_path} not found. This is just an example.")
        print("Replace with your actual h5ad file path.\n")
        return

    adata = load_h5ad(input_path)

    # Print overview of the data
    print_h5ad_overview(adata, n_examples=5)

    # Configure the filtering pipeline
    config = FilterPipelineConfig(
        pipeline_name="standard_crispr_qc",
        umi_filter=UMIFilterConfig(min_counts=15000, count_column="total_counts"),
        guide_filter=GuideFilterConfig(guide_column="pass_guide_filter", required_value=1),
        mito_filter=MitochondrialFilterConfig(
            max_pct_mt=20.0,
            mt_column="pct_counts_mt",
            enabled=True,  # Can disable individual filters
        ),
        gene_filter=GeneDetectionFilterConfig(min_genes=200, gene_column="n_genes_by_counts"),
    )

    print(f"\nPipeline configuration hash: {config.get_hash()}")

    # Run the pipeline
    adata_filtered, stats = run_pipeline(adata, config)

    # Print filtering summary
    print_filtering_summary(stats)

    # Save the filtered data
    output_path = save_filtered_h5ad(
        adata_filtered, original_path=input_path, config=config, output_dir=Path("output")
    )

    print(f"\nFiltered data saved to: {output_path}")
    print(f"Configuration saved to: {output_path.with_suffix('.json')}")


def example_batch_processing():
    """Example: Process multiple h5ad files with the same pipeline."""
    print("\n\nEXAMPLE 2: Batch Processing")
    print("=" * 80)

    # Configure pipeline (same as above, but more lenient for demo)
    config = FilterPipelineConfig(
        pipeline_name="lenient_qc",
        umi_filter=UMIFilterConfig(
            min_counts=10000,  # Lower threshold
            count_column="total_counts",
        ),
        guide_filter=GuideFilterConfig(guide_column="pass_guide_filter", required_value=1),
        # Skip mito filter by not including it
        gene_filter=GeneDetectionFilterConfig(
            min_genes=150,  # Lower threshold
            gene_column="n_genes_by_counts",
        ),
    )

    # Process all h5ad files in a directory
    input_dir = Path("data/raw")
    output_dir = Path("data/filtered")

    if not input_dir.exists():
        print(f"Note: {input_dir} not found. This is just an example.")
        print("Replace with your actual directory path.\n")
        return

    # Run batch processing
    results = batch_filter_h5ads(
        input_paths=input_dir,  # Can also pass list of specific files
        config=config,
        output_dir=output_dir,
        skip_existing=True,  # Skip files that were already processed
    )

    # Results contains detailed info about each file
    print(f"\nProcessed {len(results)} files")


def example_custom_pipeline():
    """Example: Create a custom pipeline with only specific filters."""
    print("\n\nEXAMPLE 3: Custom Pipeline")
    print("=" * 80)

    # Create a pipeline with only UMI and guide filters
    config = FilterPipelineConfig(
        pipeline_name="minimal_qc",
        umi_filter=UMIFilterConfig(min_counts=15000),
        guide_filter=GuideFilterConfig(guide_column="pass_guide_filter"),
        # mito_filter and gene_filter are None (not included)
    )

    print(f"Pipeline hash: {config.get_hash()}")
    print(f"Enabled filters: {[name for name, _ in config.get_enabled_filters()]}")

    # Save configuration for later use
    config.save_json(Path("configs/minimal_qc.json"))
    print("\nConfiguration saved to: configs/minimal_qc.json")

    # Load it back
    loaded_config = FilterPipelineConfig.load_json(Path("configs/minimal_qc.json"))
    print(f"Loaded configuration: {loaded_config.pipeline_name}")


def example_disabled_filter():
    """Example: Configure pipeline with disabled filters."""
    print("\n\nEXAMPLE 4: Selectively Disable Filters")
    print("=" * 80)

    # Create a pipeline but disable the mitochondrial filter
    config = FilterPipelineConfig(
        pipeline_name="no_mito_filter",
        umi_filter=UMIFilterConfig(min_counts=15000),
        guide_filter=GuideFilterConfig(guide_column="pass_guide_filter"),
        mito_filter=MitochondrialFilterConfig(
            max_pct_mt=20.0,
            enabled=False,  # Disabled!
        ),
        gene_filter=GeneDetectionFilterConfig(min_genes=200),
    )

    enabled = [name for name, _ in config.get_enabled_filters()]
    print(f"Enabled filters: {enabled}")
    print("Note: mito_filter is configured but disabled")


def example_different_thresholds():
    """Example: Multiple pipeline configurations with different thresholds."""
    print("\n\nEXAMPLE 5: Multiple Pipeline Configurations")
    print("=" * 80)

    # Strict QC
    strict_config = FilterPipelineConfig(
        pipeline_name="strict_qc",
        umi_filter=UMIFilterConfig(min_counts=20000),
        guide_filter=GuideFilterConfig(guide_column="pass_guide_filter"),
        mito_filter=MitochondrialFilterConfig(max_pct_mt=10.0),
        gene_filter=GeneDetectionFilterConfig(min_genes=500),
    )

    # Lenient QC
    lenient_config = FilterPipelineConfig(
        pipeline_name="lenient_qc",
        umi_filter=UMIFilterConfig(min_counts=10000),
        guide_filter=GuideFilterConfig(guide_column="pass_guide_filter"),
        mito_filter=MitochondrialFilterConfig(max_pct_mt=25.0),
        gene_filter=GeneDetectionFilterConfig(min_genes=150),
    )

    print(f"Strict config hash:  {strict_config.get_hash()}")
    print(f"Lenient config hash: {lenient_config.get_hash()}")
    print("\nThese will produce different output filenames!")

    # Both can be applied to the same input file, producing different outputs
    # data_filtered_{strict_hash}.h5ad
    # data_filtered_{lenient_hash}.h5ad


if __name__ == "__main__":
    print("filter_h5ads Package Examples")
    print("=" * 80)
    print()

    # Run examples
    # Note: Some examples require actual data files to work
    example_single_file()
    example_batch_processing()
    example_custom_pipeline()
    example_disabled_filter()
    example_different_thresholds()

    print("\n" + "=" * 80)
    print("Examples complete!")
    print("=" * 80)

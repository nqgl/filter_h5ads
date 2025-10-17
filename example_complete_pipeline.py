"""Complete filtering pipeline example with all steps.

This example demonstrates a full CRISPR screen filtering pipeline including:
1. Ensembl ID conversion (preprocessing)
2. UMI count filtering
3. Guide quality filtering
4. Mitochondrial content filtering
5. Gene detection filtering
"""

from pathlib import Path

from comlm.config.data_configs import SafeGeneSymbolRemappingRef

from filter_h5ads import (
    EnsemblConversionConfig,
    FilterPipelineConfig,
    GeneDetectionFilterConfig,
    GuideFilterConfig,
    MitochondrialFilterConfig,
    UMIFilterConfig,
    load_h5ad,
    print_filtering_summary,
    print_h5ad_overview,
    run_pipeline,
    save_filtered_h5ad,
)
from filter_h5ads.mask import mask


def main():
    """Run complete filtering pipeline on a CRISPR screen h5ad file."""
    # ==========================================================================
    # Configuration
    # ==========================================================================

    # Input/output paths
    input_path = Path("../HEK293T_filtered_dual_guide_cells.h5ad")
    output_dir = Path("data/filtered")

    # Define allowed genes mask (replace with your actual mask)
    # This should be a list of Ensembl IDs that you want to keep in the data
    # Create complete pipeline configuration
    config = FilterPipelineConfig(
        pipeline_name="complete_crispr_qc_pipeline",
        # Step 1: Ensembl ID Conversion (Preprocessing)
        # Converts gene names to Ensembl IDs and filters by allowed genes
        ensembl_conversion=EnsemblConversionConfig(
            gene_remapping_ref=SafeGeneSymbolRemappingRef.get_default_for_vcc(),
            mask_getter=lambda: mask,  # Your allowed gene set
            source_target_column="gene_target",  # Source column with gene targets
            target_column_name="target_gene",  # New column with Ensembl IDs
            enabled=True,
        ),
        # Step 2: UMI Count Filter
        # Remove low-quality cells with insufficient UMI counts
        umi_filter=UMIFilterConfig(
            min_counts=15000,  # Minimum UMI threshold (adjust based on your data)
            count_column="total_counts",
            enabled=True,
        ),
        # Step 3: Guide Quality Filter
        # Keep only cells that passed guide pairing QC
        guide_filter=GuideFilterConfig(
            guide_column="pass_guide_filter",
            required_value=1,  # Only keep cells with value = 1
            enabled=True,
        ),
        # Step 4: Mitochondrial Content Filter
        # Remove dead/dying cells with high mitochondrial percentage
        mito_filter=MitochondrialFilterConfig(
            max_pct_mt=20.0,  # Maximum 20% mitochondrial content
            mt_column="pct_counts_mt",
            enabled=True,
        ),
        # Step 5: Gene Detection Filter
        # Remove cells with too few detected genes
        gene_filter=GeneDetectionFilterConfig(
            min_genes=200,  # Minimum genes with non-zero counts
            gene_column="n_genes_by_counts",
            enabled=True,
        ),
    )

    # ==========================================================================
    # Execution
    # ==========================================================================

    print("=" * 80)
    print("COMPLETE CRISPR SCREEN FILTERING PIPELINE")
    print("=" * 80)
    print()

    # Load the data
    print(f"Loading data from: {input_path}")
    adata = load_h5ad(input_path)
    print()

    # Print overview of raw data
    print("RAW DATA OVERVIEW:")
    print_h5ad_overview(adata, n_examples=5)
    print()

    # Run the complete pipeline
    print("RUNNING FILTERING PIPELINE...")
    print(f"Pipeline: {config.pipeline_name}")
    print(f"Steps enabled: {len(config.get_enabled_filters())}")
    print()

    adata_filtered, stats = run_pipeline(adata, config)

    # Print detailed filtering summary
    print_filtering_summary(stats)

    # Save filtered data
    print("\nSAVING FILTERED DATA...")
    output_path = save_filtered_h5ad(
        adata_filtered,
        original_path=input_path,
        config=config,
        output_dir=output_dir,
    )

    print(f"\n✓ Filtered data saved to: {output_path}")
    print(f"✓ Configuration saved to: {output_path.with_suffix('.json')}")

    # Summary statistics
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"Initial cells: {adata.n_obs:,}")
    print(f"Final cells:   {adata_filtered.n_obs:,}")
    print(f"Retained:      {adata_filtered.n_obs / adata.n_obs * 100:.1f}%")
    print()
    print(f"Initial genes: {adata.n_vars:,}")
    print(f"Final genes:   {adata_filtered.n_vars:,}")
    print()
    print(f"Output: {output_path}")
    print("=" * 80)


def create_config_without_ensembl_conversion():
    """Alternative configuration without Ensembl conversion.

    Use this if your data already has Ensembl IDs or you want to skip
    the gene name conversion step.
    """
    config = FilterPipelineConfig(
        pipeline_name="standard_qc_no_conversion",
        # Skip ensembl_conversion (set to None or omit)
        ensembl_conversion=None,
        umi_filter=UMIFilterConfig(min_counts=15000),
        guide_filter=GuideFilterConfig(guide_column="pass_guide_filter"),
        mito_filter=MitochondrialFilterConfig(max_pct_mt=20.0),
        gene_filter=GeneDetectionFilterConfig(min_genes=200),
    )
    return config


def create_lenient_config():
    """Alternative configuration with more lenient thresholds.

    Use this for datasets with lower quality or when you want to retain
    more cells for downstream analysis.
    """
    mask = []  # TODO: Replace with your mask

    config = FilterPipelineConfig(
        pipeline_name="lenient_qc_pipeline",
        ensembl_conversion=EnsemblConversionConfig(
            gene_remapping_ref=SafeGeneSymbolRemappingRef.get_default_for_vcc(),
            mask_getter=lambda: mask,
            enabled=True,
        ),
        umi_filter=UMIFilterConfig(
            min_counts=10000,  # Lower threshold
            enabled=True,
        ),
        guide_filter=GuideFilterConfig(
            guide_column="pass_guide_filter",
            enabled=True,
        ),
        mito_filter=MitochondrialFilterConfig(
            max_pct_mt=25.0,  # More lenient
            enabled=True,
        ),
        gene_filter=GeneDetectionFilterConfig(
            min_genes=150,  # Lower threshold
            enabled=True,
        ),
    )
    return config


def create_strict_config():
    """Alternative configuration with stricter thresholds.

    Use this when you want higher quality cells and can afford to
    filter out more data.
    """
    mask = []  # TODO: Replace with your mask

    config = FilterPipelineConfig(
        pipeline_name="strict_qc_pipeline",
        ensembl_conversion=EnsemblConversionConfig(
            gene_remapping_ref=SafeGeneSymbolRemappingRef.get_default_for_vcc(),
            mask_getter=lambda: mask,
            enabled=True,
        ),
        umi_filter=UMIFilterConfig(
            min_counts=20000,  # Higher threshold
            enabled=True,
        ),
        guide_filter=GuideFilterConfig(
            guide_column="pass_guide_filter",
            enabled=True,
        ),
        mito_filter=MitochondrialFilterConfig(
            max_pct_mt=15.0,  # Stricter
            enabled=True,
        ),
        gene_filter=GeneDetectionFilterConfig(
            min_genes=500,  # Higher threshold
            enabled=True,
        ),
    )
    return config


def selective_filtering_example():
    """Example with selective filtering - only some steps enabled.

    This shows how to enable/disable specific filtering steps.
    """
    mask = []  # TODO: Replace with your mask

    config = FilterPipelineConfig(
        pipeline_name="selective_filtering",
        # Only run Ensembl conversion and UMI filter
        ensembl_conversion=EnsemblConversionConfig(
            gene_remapping_ref=SafeGeneSymbolRemappingRef.get_default_for_vcc(),
            mask_getter=lambda: mask,
            enabled=True,
        ),
        umi_filter=UMIFilterConfig(
            min_counts=15000,
            enabled=True,
        ),
        # Disable guide filter
        guide_filter=GuideFilterConfig(
            guide_column="pass_guide_filter",
            enabled=False,  # Disabled!
        ),
        # Skip mito and gene filters entirely (set to None)
        mito_filter=None,
        gene_filter=None,
    )
    return config


if __name__ == "__main__":
    main()

    # Uncomment to try alternative configurations:
    # config = create_config_without_ensembl_conversion()
    # config = create_lenient_config()
    # config = create_strict_config()
    # config = selective_filtering_example()

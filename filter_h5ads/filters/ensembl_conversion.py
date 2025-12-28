"""Ensembl ID conversion filter implementation.

Converts gene names to Ensembl IDs for both target genes and variable names,
and filters based on allowed gene sets for CRISPR screen data.
"""

import time
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from comlm.config.data_configs import SafeGeneSymbolRemappingRef
from loguru import logger
from pydantic import ConfigDict, Field

from filter_h5ads import mask
from filter_h5ads.filters.base import FilterStepConfig


def convert_gene_labels(
    labels: list[str],
    gene_map: dict[str, str],
    check_conflicts: bool = True,
) -> list[str]:
    """Convert gene labels to Ensembl IDs (without filtering).

    Args:
        labels: List of gene labels (mix of gene names and possibly Ensembl IDs)
        gene_map: Dictionary mapping gene names to Ensembl IDs
        check_conflicts: If True, check for conflicts between original and converted

    Returns:
        List of converted Ensembl IDs (same length as input)

    Raises:
        ValueError: If there's overlap between original and converted labels
                   (excluding "non-targeting")

    """
    original_set = set(labels)

    # Convert gene names to Ensembl IDs where possible
    converted_labels = []
    converted_set = set()

    for label in labels:
        if label in gene_map:
            ensembl_id = gene_map[label]
            converted_labels.append(ensembl_id)
            converted_set.add(ensembl_id)
        else:
            # Already an Ensembl ID or non-convertible
            converted_labels.append(label)

    # Check for conflicts (overlap between original and converted, except non-targeting)
    if check_conflicts:
        overlap = (original_set & converted_set) - {"non-targeting"}
        if overlap:
            raise ValueError(
                f"Conflict detected: Labels appear in both original and "
                f"converted sets: {overlap}. "
                "This indicates ambiguous gene identifiers that could be "
                "interpreted as both gene names and Ensembl IDs.",
            )

    return converted_labels


class EnsemblConversionConfig(FilterStepConfig):
    """Configuration for Ensembl ID conversion.

    Converts gene names to Ensembl IDs for both target genes and variable names.
    """

    mask: list[str]
    gene_remapping_ref: SafeGeneSymbolRemappingRef = Field(
        description="SafeGeneSymbolRemappingRef for gene name to Ensembl ID mapping",
    )
    source_target_column: str = Field(
        default="gene_target",
        description="Column in .obs containing gene targets to convert",
    )
    target_column_name: str = Field(
        default="target_gene",
        description="Name for new column with converted gene targets",
    )

    model_config = ConfigDict(frozen=True, extra="forbid", arbitrary_types_allowed=True)

    def mask_getter(self) -> list[str]:
        return mask.mask

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        """Apply Ensembl ID conversion to both target genes and variable names.

        This function performs two main operations:
        1. Creates a new target_gene column with converted and filtered gene targets
        2. Converts .var_names to Ensembl IDs and adds missing genes as zero columns

        Args:
            adata: Input AnnData object

        Returns:
            Tuple of (converted AnnData, statistics dictionary)

        Raises:
            ValueError: If required columns are missing or conflicts detected

        """
        logger.info("🧬 Starting Ensembl ID conversion")

        # Validate input
        if self.source_target_column not in adata.obs.columns:
            raise ValueError(
                f"Required column '{self.source_target_column}' not found in .obs. "
                f"Available columns: {list(adata.obs.columns)}",
            )

        # Open the gene remapping reference
        logger.debug("  Loading gene remapping reference...")
        t0 = time.time()
        gene_remapping = self.gene_remapping_ref.open_target()
        logger.debug(f"  ✓ Gene remapping loaded in {time.time() - t0:.2f}s")

        # Get the mask
        logger.debug("  Getting mask...")
        t0 = time.time()
        mask = self.mask_getter()
        mask_set = set(mask) if isinstance(mask, list) else mask
        logger.info(
            f"  ✓ Mask retrieved: {len(mask_set):,} genes ({time.time() - t0:.2f}s)"
        )

        # Handle backed mode: convert to memory if needed
        if adata.isbacked:
            logger.info("  ⚠️  Converting backed AnnData to memory for processing...")
            t0 = time.time()
            adata = adata.to_memory()
            logger.info(f"  ✓ Loaded into memory in {time.time() - t0:.2f}s")

        # Capture stats before modification (avoid full copy for memory efficiency)
        n_cells_initial = adata.n_obs
        n_genes_initial = adata.n_vars

        # ============================================================================
        # Step 1: Create target_gene column and filter cells
        # ============================================================================
        step1_start = time.time()
        logger.info(
            f"📍 Step 1: Converting {self.source_target_column} → "
            f"{self.target_column_name}"
        )

        # Copy and normalize the gene_target column
        logger.debug(f"  Copying {n_cells_initial:,} target genes...")
        t0 = time.time()
        target_genes = adata.obs[self.source_target_column]
        assert isinstance(target_genes, pd.Series)
        target_genes = target_genes.copy()
        logger.debug(f"  ✓ Copied in {time.time() - t0:.2f}s")

        # Normalize: "Non-Targeting" → "non-targeting"
        target_genes = target_genes.replace("Non-Targeting", "non-targeting")

        # Convert gene names to Ensembl IDs (no filtering yet)
        logger.debug(f"  Converting {len(target_genes):,} gene labels...")
        t0 = time.time()
        try:
            converted_targets = convert_gene_labels(
                labels=target_genes.tolist(),
                gene_map=gene_remapping.gene_symbol_to_ensembl_id,
                check_conflicts=True,
            )
        except ValueError as e:
            raise ValueError(f"Error converting gene targets: {e}") from e
        logger.debug(f"  ✓ Converted in {time.time() - t0:.2f}s")

        # Create new column with converted targets (same length as original)
        adata.obs[self.target_column_name] = pd.Categorical(converted_targets)

        # Filter cells: keep only those with target_gene in mask or "non-targeting"
        logger.debug("  Filtering cells by mask...")
        t0 = time.time()
        target_col = adata.obs[self.target_column_name]
        assert isinstance(target_col, pd.Series)
        valid_mask = target_col.isin(mask_set | {"non-targeting"})
        n_cells_removed_step1 = (~valid_mask).sum()

        logger.info(
            f"  Removed {n_cells_removed_step1:,} cells with targets not in mask"
        )

        adata = adata[valid_mask].copy()
        logger.debug(f"  ✓ Filtered in {time.time() - t0:.2f}s")
        logger.success(f"✅ Step 1 complete in {time.time() - step1_start:.2f}s")

        # ============================================================================
        # Step 2: Convert .var_names and add missing genes
        # ============================================================================
        step2_start = time.time()
        logger.info("📍 Step 2: Converting .var_names to Ensembl IDs")

        original_var_names = adata.var_names.tolist()
        logger.debug(f"  Processing {len(original_var_names):,} gene names...")

        # Convert variable names (no filtering yet)
        t0 = time.time()
        try:
            converted_var_names = convert_gene_labels(
                labels=original_var_names,
                gene_map=gene_remapping.gene_symbol_to_ensembl_id,
                check_conflicts=True,
            )
        except ValueError as e:
            raise ValueError(f"Error converting variable names: {e}") from e
        logger.debug(f"  ✓ Converted gene names in {time.time() - t0:.2f}s")

        # Filter to only keep genes in the mask
        # Create a boolean mask for which genes to keep
        logger.debug("  Filtering genes by mask...")
        t0 = time.time()
        genes_in_mask = [gene in mask_set for gene in converted_var_names]
        logger.debug(f"  ✓ Created gene mask in {time.time() - t0:.2f}s")

        # Subset the AnnData to only keep genes in the mask
        logger.debug(f"  Subsetting AnnData to {sum(genes_in_mask):,} genes...")
        t0 = time.time()
        adata = adata[:, genes_in_mask].copy()
        logger.debug(f"  ✓ Subsetted in {time.time() - t0:.2f}s")

        # Update var_names with the converted names
        adata.var_names = [
            converted_var_names[i] for i, keep in enumerate(genes_in_mask) if keep
        ]
        n_genes_after_conversion = len(adata.var_names)

        logger.info(
            f"  Converted .var_names: {len(original_var_names):,} → "
            f"{n_genes_after_conversion:,} genes",
        )

        # Find missing genes that need to be added
        current_genes_set = set(adata.var_names)
        missing_genes = mask_set - current_genes_set

        if missing_genes:
            logger.info(
                f"  Adding {len(missing_genes):,} missing genes as zero columns"
            )

            # Sort missing genes for deterministic ordering
            t0 = time.time()
            missing_genes_sorted = sorted(missing_genes)
            logger.debug(f"  ✓ Sorted missing genes in {time.time() - t0:.2f}s")

            # Determine if data is sparse
            is_sparse = sp.issparse(adata.X)
            logger.debug(f"  Data is {'sparse' if is_sparse else 'dense'}")

            # Create zero matrix for missing genes
            n_obs = adata.n_obs
            n_missing = len(missing_genes_sorted)

            # Get dtype, ensuring X is not None
            assert adata.X is not None, "adata.X cannot be None"
            x_dtype = adata.X.dtype

            logger.debug(f"  Creating zero matrix ({n_obs:,} × {n_missing:,})...")
            t0 = time.time()
            if is_sparse:
                # Create sparse zero matrix matching the format of existing data
                zero_matrix = sp.csr_matrix((n_obs, n_missing), dtype=x_dtype)
            else:
                # Create dense zero matrix
                zero_matrix = np.zeros((n_obs, n_missing), dtype=x_dtype)
            logger.debug(f"  ✓ Created zero matrix in {time.time() - t0:.2f}s")

            # Create AnnData for missing genes
            logger.debug("  Creating AnnData for missing genes...")
            t0 = time.time()
            new_var = pd.DataFrame(index=missing_genes_sorted)
            adata_missing = AnnData(X=zero_matrix, obs=adata.obs, var=new_var)
            logger.debug(f"  ✓ Created in {time.time() - t0:.2f}s")

            # Concatenate along gene axis using modern API
            logger.debug("  Concatenating AnnData objects...")
            t0 = time.time()
            from anndata import concat  # noqa: PLC0415

            adata = concat([adata, adata_missing], axis=1, join="outer", merge="first")
            logger.debug(f"  ✓ Concatenated in {time.time() - t0:.2f}s")

            logger.info(f"  Final gene count: {adata.n_vars:,}")

        logger.success(f"✅ Step 2 complete in {time.time() - step2_start:.2f}s")

        # Verify that mask is satisfied
        logger.debug("  Verifying mask satisfaction...")
        final_genes_set = set(adata.var_names)
        if not mask_set.issubset(final_genes_set):
            missing_after = mask_set - final_genes_set
            raise RuntimeError(
                f"Failed to satisfy mask constraint. Still missing "
                f"{len(missing_after)} genes: "
                f"{list(missing_after)[:10]}"
                f"{'...' if len(missing_after) > 10 else ''}",
            )
        logger.debug("  ✓ All mask genes present")

        # ============================================================================
        # Calculate statistics
        # ============================================================================
        n_cells_after = adata.n_obs
        n_cells_removed = n_cells_initial - n_cells_after
        pct_retained = (
            (n_cells_after / n_cells_initial * 100) if n_cells_initial > 0 else 0.0
        )

        stats = {
            "filter_name": "Ensembl Conversion",
            "n_cells_before": n_cells_initial,
            "n_cells_after": n_cells_after,
            "n_cells_removed": n_cells_removed,
            "pct_retained": round(pct_retained, 2),
            "source_column": self.source_target_column,
            "target_column": self.target_column_name,
            "n_genes_before": n_genes_initial,
            "n_genes_after": adata.n_vars,
            "n_genes_added": len(missing_genes) if missing_genes else 0,
            "n_cells_removed_step1": int(n_cells_removed_step1),
        }

        logger.success(
            f"🎉 Ensembl conversion complete: "
            f"{n_genes_initial:,} → {adata.n_vars:,} genes, "
            f"{n_cells_initial:,} → {n_cells_after:,} cells"
        )

        return adata, stats

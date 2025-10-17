"""Ensembl ID conversion filter implementation.

Converts gene names to Ensembl IDs for both target genes and variable names,
and filters based on allowed gene sets for CRISPR screen data.
"""

import logging
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from comlm.config.data_configs import SafeGeneSymbolRemappingRef
from pydantic import ConfigDict, Field

from filter_h5ads.filters.base import FilterStepConfig

logger = logging.getLogger(__name__)


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

    gene_remapping_ref: SafeGeneSymbolRemappingRef = Field(
        description="SafeGeneSymbolRemappingRef for gene name to Ensembl ID mapping",
    )
    mask_getter: Any = Field(
        description="Callable that returns list/set of allowed Ensembl IDs"
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
        logger.info("Starting Ensembl ID conversion")

        # Validate input
        if self.source_target_column not in adata.obs.columns:
            raise ValueError(
                f"Required column '{self.source_target_column}' not found in .obs. "
                f"Available columns: {list(adata.obs.columns)}",
            )

        # Open the gene remapping reference
        gene_remapping = self.gene_remapping_ref.open_target()

        # Get the mask
        mask = self.mask_getter()
        mask_set = set(mask) if isinstance(mask, list) else mask

        # Handle backed mode: convert to memory if needed
        if adata.isbacked:
            logger.info("  Converting backed AnnData to memory for processing")
            adata = adata.to_memory()

        # Capture stats before modification (avoid full copy for memory efficiency)
        n_cells_initial = adata.n_obs
        n_genes_initial = adata.n_vars

        # ============================================================================
        # Step 1: Create target_gene column and filter cells
        # ============================================================================
        logger.info(
            f"Step 1: Converting {self.source_target_column} → "
            f"{self.target_column_name}"
        )

        # Copy and normalize the gene_target column
        target_genes = adata.obs[self.source_target_column]
        assert isinstance(target_genes, pd.Series)
        target_genes = target_genes.copy()

        # Normalize: "Non-Targeting" → "non-targeting"
        target_genes = target_genes.replace("Non-Targeting", "non-targeting")

        # Convert gene names to Ensembl IDs (no filtering yet)
        try:
            converted_targets = convert_gene_labels(
                labels=target_genes.tolist(),
                gene_map=gene_remapping.gene_symbol_to_ensembl_id,
                check_conflicts=True,
            )
        except ValueError as e:
            raise ValueError(f"Error converting gene targets: {e}") from e

        # Create new column with converted targets (same length as original)
        adata.obs[self.target_column_name] = pd.Categorical(converted_targets)

        # Filter cells: keep only those with target_gene in mask or "non-targeting"
        target_col = adata.obs[self.target_column_name]
        assert isinstance(target_col, pd.Series)
        valid_mask = target_col.isin(mask_set | {"non-targeting"})
        n_cells_removed_step1 = (~valid_mask).sum()

        logger.info(
            f"  Removed {n_cells_removed_step1:,} cells with targets not in mask"
        )

        adata = adata[valid_mask].copy()

        # ============================================================================
        # Step 2: Convert .var_names and add missing genes
        # ============================================================================
        logger.info("Step 2: Converting .var_names to Ensembl IDs")

        original_var_names = adata.var_names.tolist()

        # Convert variable names (no filtering yet)
        try:
            converted_var_names = convert_gene_labels(
                labels=original_var_names,
                gene_map=gene_remapping.gene_symbol_to_ensembl_id,
                check_conflicts=True,
            )
        except ValueError as e:
            raise ValueError(f"Error converting variable names: {e}") from e

        # Filter to only keep genes in the mask
        # Create a boolean mask for which genes to keep
        genes_in_mask = [gene in mask_set for gene in converted_var_names]

        # Subset the AnnData to only keep genes in the mask
        adata = adata[:, genes_in_mask].copy()

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
            missing_genes_sorted = sorted(missing_genes)

            # Determine if data is sparse
            is_sparse = sp.issparse(adata.X)

            # Create zero matrix for missing genes
            n_obs = adata.n_obs
            n_missing = len(missing_genes_sorted)

            # Get dtype, ensuring X is not None
            assert adata.X is not None, "adata.X cannot be None"
            x_dtype = adata.X.dtype

            if is_sparse:
                # Create sparse zero matrix matching the format of existing data
                zero_matrix = sp.csr_matrix((n_obs, n_missing), dtype=x_dtype)
            else:
                # Create dense zero matrix
                zero_matrix = np.zeros((n_obs, n_missing), dtype=x_dtype)

            # Concatenate existing data with zeros
            if is_sparse:
                adata.X = sp.hstack([adata.X, zero_matrix], format="csr")
            else:
                adata.X = np.hstack([adata.X, zero_matrix])

            # Create var DataFrame for new genes
            new_var = pd.DataFrame(index=missing_genes_sorted)

            # Concatenate var DataFrames - cast to DataFrame for type checker
            var_df = adata.var
            assert isinstance(var_df, pd.DataFrame)
            adata.var = pd.concat([var_df, new_var], axis=0)

            # Update var_names to include new genes
            adata.var_names = list(adata.var.index)

            logger.info(f"  Final gene count: {adata.n_vars:,}")

        # Verify that mask is satisfied
        final_genes_set = set(adata.var_names)
        if not mask_set.issubset(final_genes_set):
            missing_after = mask_set - final_genes_set
            raise RuntimeError(
                f"Failed to satisfy mask constraint. Still missing "
                f"{len(missing_after)} genes: "
                f"{list(missing_after)[:10]}"
                f"{'...' if len(missing_after) > 10 else ''}",
            )

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

        logger.info("Ensembl conversion complete")

        return adata, stats

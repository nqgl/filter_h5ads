#!/usr/bin/env python3
"""Quick test to verify the configuration system works correctly.
Tests config creation, hashing, and serialization without requiring real data.
"""

from filter_h5ads import (
    FilterPipelineConfig,
    GeneDetectionFilterConfig,
    GuideFilterConfig,
    MitochondrialFilterConfig,
    UMIFilterConfig,
)


def test_config_creation():
    """Test creating pipeline configurations."""
    print("Testing configuration creation...")

    # Create a complete pipeline config
    config = FilterPipelineConfig(
        pipeline_name="test_pipeline",
        umi_filter=UMIFilterConfig(min_counts=15000),
        guide_filter=GuideFilterConfig(guide_column="pass_guide_filter"),
        mito_filter=MitochondrialFilterConfig(max_pct_mt=20.0),
        gene_filter=GeneDetectionFilterConfig(min_genes=200),
    )

    print(f"✓ Created config: {config.pipeline_name}")
    print(f"  Hash: {config.get_hash()}")

    # Test that config is immutable
    try:
        config.pipeline_name = "modified"  # type: ignore
        print("✗ Config should be immutable!")
    except Exception:
        print("✓ Config is immutable (as expected)")

    return config


def test_config_hash():
    """Test that identical configs produce identical hashes."""
    print("\nTesting config hashing...")

    config1 = FilterPipelineConfig(
        pipeline_name="test",
        umi_filter=UMIFilterConfig(min_counts=15000),
    )

    config2 = FilterPipelineConfig(
        pipeline_name="test",
        umi_filter=UMIFilterConfig(min_counts=15000),
    )

    config3 = FilterPipelineConfig(
        pipeline_name="test",
        umi_filter=UMIFilterConfig(min_counts=20000),  # Different threshold
    )

    hash1 = config1.get_hash()
    hash2 = config2.get_hash()
    hash3 = config3.get_hash()

    print(f"  Config 1 hash: {hash1}")
    print(f"  Config 2 hash: {hash2}")
    print(f"  Config 3 hash: {hash3}")

    if hash1 == hash2:
        print("✓ Identical configs produce identical hashes")
    else:
        print("✗ Identical configs should produce identical hashes!")

    if hash1 != hash3:
        print("✓ Different configs produce different hashes")
    else:
        print("✗ Different configs should produce different hashes!")


def test_config_serialization():
    """Test saving and loading configs."""
    print("\nTesting config serialization...")
    import tempfile
    from pathlib import Path

    config = FilterPipelineConfig(
        pipeline_name="serialize_test",
        umi_filter=UMIFilterConfig(min_counts=15000),
        guide_filter=GuideFilterConfig(guide_column="pass_guide_filter"),
    )

    # Save to temp file
    with tempfile.TemporaryDirectory() as tmpdir:
        config_path = Path(tmpdir) / "test_config.json"
        config.save_json(config_path)
        print(f"✓ Saved config to {config_path}")

        # Load it back
        loaded_config = FilterPipelineConfig.load_json(config_path)
        print(f"✓ Loaded config: {loaded_config.pipeline_name}")

        # Verify hashes match
        if config.get_hash() == loaded_config.get_hash():
            print("✓ Loaded config matches original (same hash)")
        else:
            print("✗ Loaded config should match original!")


def test_enabled_filters():
    """Test getting enabled filters."""
    print("\nTesting filter enabling/disabling...")

    config = FilterPipelineConfig(
        pipeline_name="partial_pipeline",
        umi_filter=UMIFilterConfig(min_counts=15000),
        guide_filter=GuideFilterConfig(guide_column="pass_guide_filter", enabled=False),
        mito_filter=None,  # Not included at all
        gene_filter=GeneDetectionFilterConfig(min_genes=200),
    )

    enabled = config.get_enabled_filters()
    enabled_names = [name for name, _ in enabled]

    print(f"  Enabled filters: {enabled_names}")

    if "umi_filter" in enabled_names and "gene_filter" in enabled_names:
        print("✓ Enabled filters detected correctly")
    else:
        print("✗ Should have umi_filter and gene_filter enabled")

    if "guide_filter" not in enabled_names and "mito_filter" not in enabled_names:
        print("✓ Disabled/missing filters excluded correctly")
    else:
        print("✗ guide_filter and mito_filter should not be enabled")


def test_validation():
    """Test config validation."""
    print("\nTesting config validation...")

    # Test invalid values
    try:
        _ = UMIFilterConfig(min_counts=-100)  # Negative not allowed
        print("✗ Should not allow negative min_counts")
    except Exception:
        print("✓ Rejects negative min_counts")

    try:
        _ = MitochondrialFilterConfig(max_pct_mt=150)  # > 100% not allowed
        print("✗ Should not allow max_pct_mt > 100")
    except Exception:
        print("✓ Rejects max_pct_mt > 100")

    # Test valid edge cases
    try:
        _ = UMIFilterConfig(min_counts=0)  # Zero is valid
        print("✓ Allows min_counts=0")
    except Exception:
        print("✗ Should allow min_counts=0")


if __name__ == "__main__":
    print("=" * 80)
    print("filter_h5ads Configuration Tests")
    print("=" * 80)
    print()

    test_config_creation()
    test_config_hash()
    test_config_serialization()
    test_enabled_filters()
    test_validation()

    print("\n" + "=" * 80)
    print("All configuration tests complete!")
    print("=" * 80)

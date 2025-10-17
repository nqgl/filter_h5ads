#!/bin/bash
# Example CLI usage for filter_h5ads

echo "filter_h5ads CLI Examples"
echo "========================="
echo ""

# Example 1: Inspect a file
echo "1. Inspect an h5ad file:"
echo "   filter-h5ads inspect data.h5ad"
echo ""

# Example 2: Inspect with more examples
echo "2. Inspect with more examples:"
echo "   filter-h5ads inspect data.h5ad --n-examples 10"
echo ""

# Example 3: Get help
echo "3. Get help:"
echo "   filter-h5ads --help"
echo "   filter-h5ads inspect --help"
echo ""

# Example 4: Filter a single file
echo "4. Filter a single file using a config:"
echo "   filter-h5ads filter data.h5ad config.json"
echo ""

# Example 5: Filter with custom output directory
echo "5. Filter with custom output directory:"
echo "   filter-h5ads filter data.h5ad config.json --output-dir filtered/"
echo ""

# Example 6: Batch processing
echo "6. Batch process all h5ad files in a directory:"
echo "   filter-h5ads batch data/raw/ config.json --output-dir data/filtered/"
echo ""

# Example 7: Batch processing without skipping existing
echo "7. Batch process, reprocessing existing files:"
echo "   filter-h5ads batch data/raw/ config.json --output-dir data/filtered/ --no-skip-existing"
echo ""

echo "========================="
echo "Note: Replace data.h5ad and config.json with your actual file paths"

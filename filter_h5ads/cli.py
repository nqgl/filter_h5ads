"""
Command-line interface for filter_h5ads package.
"""

import sys
from pathlib import Path

import click

from filter_h5ads.inspection import print_h5ad_overview
from filter_h5ads.io import load_h5ad


@click.group()
@click.version_option(version="0.1.0")
def cli():
    """filter_h5ads: Filter and inspect h5ad files from CRISPR screens."""
    pass


@cli.command()
@click.argument("h5ad_path", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "--n-examples",
    "-n",
    default=5,
    type=int,
    help="Number of unique examples to show for each column",
)
@click.option(
    "--backed",
    is_flag=True,
    help="Load h5ad file in backed mode (memory efficient)",
)
def inspect(h5ad_path: Path, n_examples: int, backed: bool):
    """Inspect an h5ad file and print comprehensive overview.

    \b
    Example:
        filter-h5ads inspect data.h5ad
        filter-h5ads inspect data.h5ad --n-examples 10
        filter-h5ads inspect data.h5ad --backed
    """
    try:
        adata = load_h5ad(h5ad_path, backed=backed)
        print_h5ad_overview(adata, n_examples=n_examples)
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.argument("h5ad_path", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("config_path", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(file_okay=False, path_type=Path),
    help="Output directory (defaults to same as input)",
)
@click.option(
    "--backed",
    is_flag=True,
    help="Load h5ad file in backed mode (memory efficient)",
)
def filter_h5ad(h5ad_path: Path, config_path: Path, output_dir: Path | None, backed: bool):
    """Filter an h5ad file using a configuration file.

    \b
    Example:
        filter-h5ads filter data.h5ad config.json
        filter-h5ads filter data.h5ad config.json --output-dir filtered/
        filter-h5ads filter data.h5ad config.json --backed
    """
    from filter_h5ads.config import FilterPipelineConfig
    from filter_h5ads.inspection import print_filtering_summary
    from filter_h5ads.io import save_filtered_h5ad
    from filter_h5ads.pipeline import run_pipeline, validate_pipeline_config

    try:
        # Load data
        click.echo(f"Loading {h5ad_path}...")
        adata = load_h5ad(h5ad_path, backed=backed)

        # Load config
        click.echo(f"Loading configuration from {config_path}...")
        config = FilterPipelineConfig.load_json(config_path)

        # Validate config
        errors = validate_pipeline_config(config, adata)
        if errors:
            click.echo("Configuration validation failed:", err=True)
            for error in errors:
                click.echo(f"  • {error}", err=True)
            sys.exit(1)

        # Run pipeline
        click.echo(f"Running pipeline: {config.pipeline_name}")
        adata_filtered, stats = run_pipeline(adata, config)

        # Print summary
        print_filtering_summary(stats)

        # Save results
        output_path = save_filtered_h5ad(
            adata_filtered, original_path=h5ad_path, config=config, output_dir=output_dir
        )

        click.echo(f"\n✓ Filtered data saved to: {output_path}")
        click.echo(f"✓ Configuration saved to: {output_path.with_suffix('.json')}")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.argument("input_dir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.argument("config_path", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "--output-dir",
    "-o",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="Output directory for filtered files",
)
@click.option(
    "--skip-existing/--no-skip-existing",
    default=True,
    help="Skip files that have already been processed",
)
@click.option(
    "--backed",
    is_flag=True,
    help="Load h5ad files in backed mode (memory efficient)",
)
def batch(input_dir: Path, config_path: Path, output_dir: Path, skip_existing: bool, backed: bool):
    """Batch filter multiple h5ad files in a directory.

    \b
    Example:
        filter-h5ads batch data/raw/ config.json --output-dir data/filtered/
        filter-h5ads batch data/raw/ config.json --output-dir data/filtered/ --backed
    """
    from filter_h5ads.config import FilterPipelineConfig
    from filter_h5ads.io import batch_filter_h5ads

    try:
        # Load config
        click.echo(f"Loading configuration from {config_path}...")
        config = FilterPipelineConfig.load_json(config_path)

        # Run batch processing
        click.echo(f"Processing h5ad files in {input_dir}...")
        batch_filter_h5ads(
            input_paths=input_dir,
            config=config,
            output_dir=output_dir,
            skip_existing=skip_existing,
            backed=backed,
        )

        # Summary is printed by batch_filter_h5ads
        click.echo("\n✓ Batch processing complete!")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    cli()

"""Obs column transform step.

Creates a new `.obs` column by applying a sequence of safe, serializable operations
to one or more input columns. This is intended for light metadata normalization
and feature derivation without using `eval()`.
"""

from __future__ import annotations

import ast
from collections.abc import Iterable, Iterator
from typing import Annotated, Any, Literal, Union

import pandas as pd
from anndata import AnnData
from loguru import logger
from pydantic import BaseModel, ConfigDict, Field

from filter_h5ads.filters.base import FilterStepConfig
from filter_h5ads.utils import calculate_filter_stats, validate_h5ad_columns


class _ValueOp(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    op: str

    def run(self, value: Any) -> Any:  # pragma: no cover - exercised via subclasses
        raise NotImplementedError


class ParsePythonLiteralOp(_ValueOp):
    """Parse a Python literal using `ast.literal_eval` (safe subset of Python syntax)."""

    op: Literal["parse_python_literal"] = "parse_python_literal"
    only_if_str: bool = Field(
        default=True,
        description="If True, only parses when the value is a string/bytes.",
    )

    def run(self, value: Any) -> Any:
        if self.only_if_str and not isinstance(value, (str, bytes)):
            return value
        if isinstance(value, bytes):
            value = value.decode("utf-8")
        if not isinstance(value, str):
            raise TypeError(
                f"parse_python_literal expected str/bytes, got {type(value).__name__}"
            )
        return ast.literal_eval(value)


class FirstOp(_ValueOp):
    """Take the first element of a list/tuple."""

    op: Literal["first"] = "first"

    def run(self, value: Any) -> Any:
        if not isinstance(value, (list, tuple)):
            raise TypeError(f"first expected list/tuple, got {type(value).__name__}")
        if len(value) == 0:
            raise ValueError("first expected a non-empty list/tuple")
        return value[0]


class IndexOp(_ValueOp):
    """Index into a list/tuple."""

    op: Literal["index"] = "index"
    index: int = Field(description="Index to select")

    def run(self, value: Any) -> Any:
        if not isinstance(value, (list, tuple)):
            raise TypeError(f"index expected list/tuple, got {type(value).__name__}")
        return value[self.index]


class PickIndicesOp(_ValueOp):
    """Pick multiple indices from a list/tuple and return them as a tuple."""

    op: Literal["pick_indices"] = "pick_indices"
    indices: list[int] = Field(min_length=1, description="Indices to pick, in order")

    def run(self, value: Any) -> Any:
        if not isinstance(value, (list, tuple)):
            raise TypeError(
                f"pick_indices expected list/tuple, got {type(value).__name__}"
            )
        return tuple(value[i] for i in self.indices)


class JoinOp(_ValueOp):
    """Join elements of a list/tuple into a string."""

    op: Literal["join"] = "join"
    sep: str = Field(default=" ", description="Separator used to join elements")

    def run(self, value: Any) -> Any:
        if not isinstance(value, (list, tuple)):
            raise TypeError(f"join expected list/tuple, got {type(value).__name__}")
        return self.sep.join("" if v is None else str(v) for v in value)


class ToStringOp(_ValueOp):
    """Convert the value to a string."""

    op: Literal["to_string"] = "to_string"

    def run(self, value: Any) -> Any:
        return "" if value is None else str(value)


ValueOp = Annotated[
    Union[
        ParsePythonLiteralOp,
        FirstOp,
        IndexOp,
        PickIndicesOp,
        JoinOp,
        ToStringOp,
    ],
    Field(discriminator="op"),
]


class ObsColumnTransformConfig(FilterStepConfig):
    """Derive a new `.obs` column by transforming one or more existing columns."""

    input_columns: list[str] = Field(
        min_length=1,
        description="One or more `.obs` columns to read as the input value.",
    )
    input_mode: Literal["value", "tuple", "dict"] = Field(
        default="value",
        description=(
            "How to construct the initial value from input_columns. "
            "'value' requires exactly one input column; "
            "'tuple' produces a tuple of input column values; "
            "'dict' produces a dict mapping column name to value."
        ),
    )
    operations: list[ValueOp] = Field(
        min_length=1,
        description="Operations applied sequentially to derive the output value.",
    )
    output_column: str = Field(description="Name of the new `.obs` column to create.")
    overwrite: bool = Field(
        default=False,
        description="If True, allow overwriting an existing output_column.",
    )
    on_error: Literal["raise", "set_null"] = Field(
        default="raise",
        description="Behavior when an operation fails for a row.",
    )
    output_dtype: Literal["string", "category", "object"] = Field(
        default="string",
        description="Pandas dtype to use for the output column.",
    )

    model_config = ConfigDict(frozen=True, extra="forbid")

    def _apply_ops(self, value: Any) -> Any:
        current = value
        for op in self.operations:
            current = op.run(current)
        return current

    def apply(self, adata: AnnData) -> tuple[AnnData, dict[str, Any]]:
        """Create a new `.obs` column with derived values."""
        logger.info(
            "Applying obs column transform: "
            f"{self.input_columns} → {self.output_column}"
        )

        validate_h5ad_columns(adata, self.input_columns, location="obs")
        if getattr(adata, "is_view", False):
            adata = adata.copy()

        if self.output_column in adata.obs.columns and not self.overwrite:
            raise ValueError(
                f"Output column '{self.output_column}' already exists in .obs. "
                "Set overwrite=True to replace it."
            )

        if self.input_mode == "value" and len(self.input_columns) != 1:
            raise ValueError("input_mode='value' requires exactly one input column")

        df = adata.obs
        if self.input_mode == "value":
            series = df[self.input_columns[0]]
            assert isinstance(series, pd.Series)
            base_values_iter: Iterable[Any] = series.array
        elif self.input_mode == "tuple":
            cols = [df[c].array for c in self.input_columns]
            base_values_iter = zip(*cols, strict=True)
        else:  # dict
            cols = {c: df[c].array for c in self.input_columns}
            n = len(df)

            def _dict_iter() -> Iterator[dict[str, Any]]:
                for i in range(n):
                    yield {c: cols[c][i] for c in self.input_columns}

            base_values_iter = _dict_iter()

        outputs: list[Any] = []
        n_errors = 0
        for base_value in base_values_iter:
            try:
                outputs.append(self._apply_ops(base_value))
            except Exception:
                if self.on_error == "raise":
                    raise
                n_errors += 1
                outputs.append(None)

        out_series = pd.Series(outputs, index=df.index)
        if self.output_dtype == "string":
            out_series = out_series.astype("string")
        elif self.output_dtype == "category":
            out_series = out_series.astype("category")

        adata.obs[self.output_column] = out_series

        stats = calculate_filter_stats(
            adata_before=adata,
            adata_after=adata,
            filter_name="Obs Column Transform",
        )
        stats["input_columns"] = list(self.input_columns)
        stats["input_mode"] = self.input_mode
        stats["output_column"] = self.output_column
        stats["overwrite"] = self.overwrite
        stats["output_dtype"] = self.output_dtype
        stats["n_row_errors"] = int(n_errors)

        logger.info(
            f"Obs column transform complete: wrote '{self.output_column}' "
            f"(row errors: {n_errors})"
        )
        return adata, stats

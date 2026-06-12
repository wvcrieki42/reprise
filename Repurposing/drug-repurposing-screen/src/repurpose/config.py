"""Typed configuration loader."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import yaml


@dataclass
class Config:
    raw: dict[str, Any]
    root: Path

    # convenience accessors -------------------------------------------------
    @property
    def mode(self) -> str:
        return self.raw.get("mode", "demo")

    def path(self, key: str) -> Path:
        p = Path(self.raw["paths"][key])
        return p if p.is_absolute() else (self.root / p)

    def get(self, *keys: str, default: Any = None) -> Any:
        node: Any = self.raw
        for k in keys:
            if not isinstance(node, dict) or k not in node:
                return default
            node = node[k]
        return node

    @property
    def outdir(self) -> Path:
        return self.path("outdir")


def load_config(path: str | Path) -> Config:
    path = Path(path).resolve()
    with open(path) as fh:
        raw = yaml.safe_load(fh)
    cfg = Config(raw=raw, root=path.parent)
    cfg.outdir.mkdir(parents=True, exist_ok=True)
    return cfg

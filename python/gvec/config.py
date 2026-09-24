# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gvec.core.config - contains the configuration classes for the GVEC library"""

# === Imports === #

from abc import ABC
from collections.abc import Mapping
import warnings

# === Config class for picard_current === #


class ConfigABC(ABC):
    aliases = {}

    @classmethod
    def from_dict(cls, config: Mapping | None):
        if config is None:
            return cls()

        for key in config:
            if key in cls.aliases:
                warnings.warn(
                    f"Config key {key} is deprecated. Use {cls.aliases[key]} instead.",
                    DeprecationWarning,
                )
        config = {cls.aliases.get(k, k): v for k, v in config.items()}

        return cls(**config)

    def __getitem__(self, key):
        if key in self.keys():
            return getattr(self, key)
        raise KeyError(key)

    @classmethod
    def default(cls, key):
        import inspect

        return inspect.signature(cls.__init__).parameters[key].default

    @classmethod
    def keys(cls):
        return cls.__init__.__annotations__.keys()

    def to_dict(self):
        return {key: self[key] for key in self.keys()}

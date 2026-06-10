import yaml
from types import SimpleNamespace

with open("config.yaml") as f:
    params = SimpleNamespace(**yaml.safe_load(f))
    
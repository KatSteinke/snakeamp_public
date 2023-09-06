"""Load default pipeline config."""

import pathlib

import yaml

# import check_config - TODO: what config checks do we need in this one?

default_config_file = pathlib.Path(__file__).parent / "config" / "pipeline_routine.yaml"

with open(default_config_file, "r", encoding="utf-8") as default_config:
    WORKFLOW_DEFAULT_CONF = yaml.safe_load(default_config)

#check_config.check_config(WORKFLOW_DEFAULT_CONF)

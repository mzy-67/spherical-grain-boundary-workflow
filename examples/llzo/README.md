# LLZO reproduction example

This example preserves the chemistry, atom-type order, origin anchor, transport carrier, and temperatures used by the original LLZO study. Update the potential path and cluster settings before running it.

```bash
cp examples/llzo/config.yaml configs/llzo.local.yaml
python scripts/submit_gb_workflow.py --config configs/llzo.local.yaml --start-row 0 --stop-row 1 --dry-run
```

# Single-GB smoke test

This one-row LLZO example exercises the same generic configuration schema as a production run. Set a valid potential and cluster configuration, then run:

```bash
python scripts/submit_gb_workflow.py --config examples/single_gb/config.yaml --start-row 0 --stop-row 1 --dry-run
```

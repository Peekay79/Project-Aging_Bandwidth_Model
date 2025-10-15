Reproduction commands:

1) Install dependencies
   pip install -r requirements.txt

2) Download data (idempotent)
   python scripts/download_data.py  # set env vars as needed

3) Prepare data
   python scripts/prepare_data.py

4) Run main pipeline
   python scripts/run_pipeline.py

5) Inflection analysis
   python scripts/inflection.py

6) Network overlay hotspots
   python scripts/network_overlay.py  # hotspots run from pipeline; script available for utilities

# Curvature Demo — run instructions

Open in Google Colab (replace USERNAME/REPO with your GitHub path):

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/DavidEngland/ABL/blob/main/notebooks/Curvature_Demo.ipynb)

Notes
- The notebook includes a top cell that installs minimal requirements when run in Colab.
- If running locally: `pip install -r ../requirements.txt` then `jupyter lab`.
- The notebook writes a tiny `sample_profiles.nc` for demonstrations; delete if you prefer not to keep generated files.

Quick troubleshooting
- If matplotlib style `seaborn-*` fails, install seaborn (`pip install seaborn`) or change the style to a matplotlib default (`'ggplot'`).
- If ipywidgets controls don't render, enable widgets extension in JupyterLab: `jupyter labextension install @jupyter-widgets/jupyterlab-manager` (or use classic notebook).

Suggested collaboration flow
1. Open notebook in Colab for interactive exploration.
2. Run "Core functions" cell then "Interactive controls" cell.
3. Use synthetic NetCDF cell to generate small test data, or point to your local NetCDFs.

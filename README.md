# LIFT
This repository contains the implementation and experimental code for **LIFT**, a latent flow matching framework for probabilistic cell fate prediction from single-cell snapshots and lineage-tracing information.

LIFT combines a conditional latent-variable model with a Mean Flow decoder to represent one-to-many cell-state transitions and enable efficient one-step conditional generation.

## Repository Structure

```
LIFT/
│
├── data/
│   └── (Dataset files should be downloaded manually)
│
├── toy/
│   └── toy_lift.ipynb
│   └── toy_flow_matching.ipynb
│   └── toy_meanflow.ipynb
│   └── toy_score_matching.ipynb
│
├── real/
│   └── real_LIFT.ipynb
│   └── real_CFM.ipynb
│   └── real_MF.ipynb
│
├── utils.py
│
└── README.md
```

- `data/`: Directory for storing datasets used in the experiments. The original datasets are not included in this repository and must be downloaded separately.
- `toy/`: Experiments on two-dimensional synthetic branching datasets, including model training, conditional generation, and evaluation.
- `real/`: Experiments on real-world lineage-traced single-cell datasets, including mouse hematopoiesis and fibroblast reprogramming.
- `utils.py`: Shared utility functions used across the experiments, including common computational functions.

## Data Sources

The real-data experiments use two publicly available lineage-tracing datasets. **The** `data/` **directory does not contain the original datasets. Users must download the required data from the original sources and place the downloaded files in this directory before running the corresponding experiments.**

### 1. Mouse Hematopoiesis: LARRY Dataset

The hematopoiesis experiment uses the LARRY lineage-tracing dataset from:

Weinreb, C., Rodriguez-Fraticelli, A., Camargo, F. D., and Klein, A. M. (2020). *Lineage tracing on transcriptional landscapes links state to fate during differentiation*. Science, 367, eaaw3381.

- **Original paper:** https://doi.org/10.1126/science.aaw3381
- **GEO accession:** [GSE140802](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140802)
- **Processed data and supplementary resources:** [AllonKleinLab/paper-data](https://github.com/AllonKleinLab/paper-data/tree/master/Lineage_tracing_on_transcriptional_landscapes_links_state_to_fate_during_differentiation)

The experiment in this repository considers the transition from Day 2 to Day 6 and predicts monocyte (Mo) and neutrophil (Neu) cell fates.

### 2. Fibroblast Reprogramming: CellTagging Dataset

The reprogramming experiment uses the CellTagging dataset from:

Biddy, B. A., Kong, W., Kamimoto, K., Guo, C., Waye, S. E., Sun, T., and Morris, S. A. (2018). *Single-cell mapping of lineage and identity in direct reprogramming*. Nature, 564, 219–224.

- **Original paper:** https://doi.org/10.1038/s41586-018-0744-4
- **GEO accession:** [GSE99915](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99915)
- **CellTagging resources:** https://morrislab.io/celltagging-resources/
- **Original data-processing code:** https://github.com/morris-lab/CellTagWorkflow

The experiment in this repository considers the transition from Day 15 to Day 28 and predicts induced endoderm progenitor (iEP) and dead-end cell fates.

## Quick Start: Synthetic Experiments

The easiest way to get started with LIFT is to run the two-dimensional synthetic experiments.

1. Clone this repository:

   ```
   git clone <REPOSITORY_URL>
   cd LIFT
   ```

2. Prepare a Python environment.

3. Open the following notebook:

   ```
   toy/toy_lift.ipynb
   ```

4. Run all cells in `toy_lift.ipynb` sequentially.

This notebook provides a self-contained entry point for reproducing the synthetic branching experiments, including model training, conditional sampling, and evaluation of population reconstruction and branching probabilities.

No external biological datasets are required for the toy experiments.

## Real-Data Experiments

The experiments on real-world datasets are organized in the `real/` directory.

Before running these experiments:

1. Download the LARRY and/or CellTagging datasets from the original sources listed above.
2. Place the required data files in the `data/` directory.
3. Check the corresponding scripts and notebooks in `real/` for data-loading paths and experiment-specific settings.
4. Run the relevant experiments to reproduce the hematopoiesis and fibroblast reprogramming results.

## Citation

If you find this repository useful in your research, please cite the corresponding LIFT paper:

```
@article{yang2026lift,
  title={LIFT: Latent Flow Matching for Probabilistic Cell Fate Dynamics with Single Cell Snapshots},
  author={Bowen,Yang and Junda, Ying and Yuxuan, Wang and Peijie, Zhou and Lei, Zhang},
  year={2026}
}
```

## License

## License

This project is licensed under the MIT License.
See the [LICENSE](LICENSE) file for details

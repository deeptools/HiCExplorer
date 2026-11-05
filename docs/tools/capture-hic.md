# Capture Hi-C

| Tool | Description |
|---|---|
| [chicQualityControl](chicQualityControl.md) | Quality control for capture Hi-C viewpoints: checks sparsity and removes viewpoints that are too sparse. |
| [chicViewpointBackgroundModel](chicViewpointBackgroundModel.md) | Computes the background model used by the other capture Hi-C tools. |
| [chicViewpoint](chicViewpoint.md) | Computes one viewpoint file per sample per reference point, based on the background model. |
| [chicSignificantInteractions](chicSignificantInteractions.md) | Detects significant interactions per viewpoint from the background model. |
| [chicAggregateStatistic](chicAggregateStatistic.md) | Aggregates viewpoint data of two samples into targets for the differential test. |
| [chicDifferentialTest](chicDifferentialTest.md) | Tests for differential interactions between two samples (chi2 or Fisher). |
| [chicPlotViewpoint](chicPlotViewpoint.md) | Plots a viewpoint together with its background model, significant and differential regions. |
| [chicExportData](chicExportData.md) | Exports the data stored in the capture Hi-C intermediate HDF5 files to text or bigWig. |
| [chicChicagoBackgroundModel](chicChicagoBackgroundModel.md) | Estimates the CHiCAGO background model from `.chinput` files or matrices. |
| [chicChicagoScores](chicChicagoScores.md) | Computes CHiCAGO p-values and weighted scores per interaction. |
| [chicChicagoSignificantInteractions](chicChicagoSignificantInteractions.md) | Selects the interactions above the CHiCAGO score threshold. |
| [chicChicagoPlotViewpoint](chicChicagoPlotViewpoint.md) | Plots CHiCAGO viewpoints as scatter or arcs, or writes a pyGenomeTracks links file. |

hicAggregateContacts --matrix ../Li_et_al_2015.cool --BED test_regions.bed \
    --outFileName master_aggregate.png \
    --numberOfBins 30 --range 50000:900000

hicAggregateContacts --matrix ../Li_et_al_2015.cool --BED test_regions.bed \
    --outFileName master_aggregate_hclust4.png \
    --numberOfBins 30 --range 50000:900000 --hclust 4 \
    --diagnosticHeatmapFile master_heatmap.png --clusterOnDiagonal

hicAggregateContacts --matrix ../Li_et_al_2015.cool --BED test_regions.bed \
    --outFileName master_aggregate_3d.png \
    --numberOfBins 30 --range 50000:900000 --hclust 2 \
    --plotType 3d \
    --BED2 test_regions.bed

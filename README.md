# EPPE pressure sensor placement experiments

MATLAB R2024a. Open this folder in MATLAB. The four fixed instances contain the hydraulic coverage matrices, installation constraints, common starting layouts and candidate priorities used in the article.

```matlab
main(1)                    % six algorithms, one paired run per case
main                       % six algorithms, 20 paired runs per case
ablation                   % component experiment and random-exchange control
analyze_results('archive')  % reported statistics from retained records
plot_results('archive')     % Figures 1-4 from retained records
diagnostics('verify')      % archived bounds and independent outage checks
diagnostics('budget',1)    % one paired run per reduced-budget case
diagnostics('threshold',1) % one paired run per changed-threshold case
diagnostics('sensitivity') % reevaluate fixed layouts under changed conditions
```

`main` and `ablation` reuse completed files in `results/runs`. New runs and their analysis are saved separately from the retained article records. Use `plot_results('new')` after completing the main and component runs. `EPPE` in `algorithms` is the shared PPE search, with switches for guided addition, protected removal and random exchange. Each algorithm has 4000 complete objective queries, including cache lookups; proxy calculations are counted separately.

`diagnostics` provides the nominal-floor model, budget and threshold analyses. Time-limited MILP reruns can yield different bounds; they do not replace archived Table 6 values. Budget/threshold archived adjusted probabilities retain the original prespecified correction families. Newly executed reduced comparisons must not be substituted for those probabilities.

Hydraulic optimization uses the supplied matrices and does not need a running EPANET installation. For a hydraulic rebuild, install the official EPANET MATLAB toolkit and call `prepare_data(toolkitFolder)`. This writes new matrices under `results/hydraulics` without changing the fixed study cases. Optimization Toolbox is needed for model/budget diagnostics; Statistics and Machine Learning Toolbox is used for quantiles and new spectral partitions. Figures are saved as 600 dpi RGB TIFF and MATLAB FIG files.

The comparator implementations retain their algorithmic source comments. EROA is the common-budget adaptation of the supplied IROA2/Pos_init implementation; SOA uses Dhiman and Kumar (2019), equations (5)-(14). The shared warm start and query cap are study-specific adaptations. Network and algorithm source ownership is unchanged.

Validation: 36 complete search replays matched retained layouts, objective values, 4000-point histories, score/cache counts and proxy counts. Main and component bootstrap statistics matched the retained records. Main, ablation, figure, budget, threshold and outage-check entry points were executed. The optional full hydraulic rebuild was not rerun as part of package cleanup.

Modena model source: WaterFutures/WaterBenchmarkHub, https://waterfutures.github.io/WaterBenchmarkHub/benchmarks/network-Modena.html. Model SHA-256: 65AD9991753C9280A5871AA247ED178BE18713F93136ED88F929330FEEB3AB32. The distributed INP files are retained without geometry or hydraulic edits.

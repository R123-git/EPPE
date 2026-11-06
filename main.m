function records=main(runs)
% Six-algorithm comparison. main(1) runs one paired replicate per case.
c=config();if nargin<1,runs=c.runs;end
records=run_experiments(c.main,runs,c);
analyze_results('new');
end

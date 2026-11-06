function records=ablation(runs)
% Four component variants, plus the matched random-exchange control.
c=config();if nargin<1,runs=c.runs;end
records=run_experiments(c.ablation,runs,c);
analyze_results('new');
end

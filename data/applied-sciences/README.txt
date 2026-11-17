File S1. Supporting records for the computational-cost and optimization-bound results

These records accompany Incremental Pressure Sensor Placement in Water Distribution Networks: Trade-Offs between Nominal and Single-Outage Coverage.

Table 6
  table6_timing_replays.csv: Five timing replays for each displayed method in each scenario. Seconds is the search time for the replay. FE counts all objective queries, including cache hits. ActualScores and CacheHits split that total. ProxyCalls counts proxy visits.
  table6_timing_summary.csv: Archived summary for the displayed methods. ProxyCalls is the mean per run. MedianSeconds is the median runtime. The archived method code PPE_M2 is displayed as R1-only in both files.

Table 7
  table7_bounds.csv: Full-precision nominal optima, feasible-layout values, certified statuses and bounds. L0 and L2 permit nominal coverage losses of approximately zero and at most two percentage points, respectively. The F and Q fields here are percentages. The RelaxationGain fields are percentage points. Certified equals 1 when the solver certified the corresponding optimum and 0 otherwise.
  table7_solver_calls.csv: Initial and extension MILP calls, with the solver time limits and elapsed times in seconds. Exitflag records the MATLAB intlinprog status. The final comparison is reported in table7_bounds.csv.
  table7_single_outage_checks.csv: Independently evaluated remaining coverages for each single-sensor outage of the final L0 and L2 feasible layouts. TypeL0orL2 code 1 denotes L0 and 2 denotes L2. The FailedSensor column lists the installed sensor IDs one by one, and Existing = 1 denotes a retained instrument. Fpercent is the post-outage weighted coverage, CNpercent the node coverage and CWpercent the demand coverage. The minimum Fpercent in each scenario and layout gives its Q.

The MATLAB code, fixed input files and other archived records are available at https://github.com/R123-git/EPPE (revision 871f3b2).

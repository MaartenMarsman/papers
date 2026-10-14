#!/bin/zsh
# Rerun of the sweeps that depend on the G-Wishart-envelope rejection samplers,
# after the completion-weight fix of 2026-10-13
# (notes/rejection-sampler-completion-weight.md). Runs the named scripts in
# sequence (each uses 12 cores); default: all three. Launch detached, e.g.
#   nohup zsh scripts/rerun_realized_marginals_2026-10-13.sh > R/results/rerun_2026-10-13_driver.log 2>&1 &
# History: partial_correlation_sweep started 16:44 under the sequential driver;
# at 19:45 the driver was stopped and marginal_offdiag_sweep plus
# rho_uniform_calibration were started alongside it (part-2 driver log).
ROOT="${SPIKESLAB_ROOT:-$(pwd)}"
cd "$ROOT" || exit 1
export R_PROFILE_USER=/dev/null
unset SPIKESLAB_QS
scripts=("$@")
[[ ${#scripts} -eq 0 ]] && scripts=(partial_correlation_sweep marginal_offdiag_sweep rho_uniform_calibration)
for s in "${scripts[@]}"; do
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] start $s"
  Rscript "R/scripts/$s.R" > "R/results/$s.log" 2>&1
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] end $s (exit $?)"
done
echo "[$(date '+%Y-%m-%d %H:%M:%S')] all done"

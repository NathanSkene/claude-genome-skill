#!/usr/bin/env bash
set -euo pipefail

# Monitor nf-core/raredisease pipeline on Imperial HPC
# Usage: ./monitor.sh [queue|log|progress|disk|all]

OUTDIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/results"
WORKDIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/nf-work"
REFDIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/references"
FASTQDIR="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/GenomeAnalysis/fastq"

CMD="${1:-all}"

show_queue() {
    echo "=== PBS Queue Status ==="
    echo ""
    echo "--- Your jobs ---"
    qstat -u ${HPC_USER:-YOUR_USERNAME} 2>/dev/null || \
        /opt/pbs/bin/qstat -u ${HPC_USER:-YOUR_USERNAME} 2>/dev/null || \
        echo "qstat not available (not on HPC login node?)"
    echo ""
    echo "--- Job counts by status ---"
    qstat -u ${HPC_USER:-YOUR_USERNAME} -f 2>/dev/null | grep "job_state" | sort | uniq -c | sort -rn || true
    echo ""
}

show_log() {
    echo "=== Nextflow Log (last 50 lines) ==="
    echo ""
    if [[ -f .nextflow.log ]]; then
        tail -50 .nextflow.log
    elif [[ -f "${OUTDIR}/pipeline_info/nextflow_stdout.log" ]]; then
        tail -50 "${OUTDIR}/pipeline_info/nextflow_stdout.log"
    else
        echo "No Nextflow log found."
        echo "Expected locations:"
        echo "  .nextflow.log (current directory)"
        echo "  ${OUTDIR}/pipeline_info/nextflow_stdout.log"
    fi
    echo ""
}

show_progress() {
    echo "=== Pipeline Progress ==="
    echo ""

    # Check trace file for process completion stats
    TRACE="${OUTDIR}/pipeline_info/trace.txt"
    if [[ -f "${TRACE}" ]]; then
        echo "--- Process completion summary ---"
        tail -n +2 "${TRACE}" | \
            awk -F'\t' '{status[$6]++; proc[$4][$6]++} END {
                printf "%-40s %8s %8s %8s\n", "PROCESS", "COMPLETED", "FAILED", "CACHED";
                printf "%-40s %8s %8s %8s\n", "-------", "---------", "------", "------";
                for (p in proc) {
                    printf "%-40s %8d %8d %8d\n", p,
                        (proc[p]["COMPLETED"] ? proc[p]["COMPLETED"] : 0),
                        (proc[p]["FAILED"] ? proc[p]["FAILED"] : 0),
                        (proc[p]["CACHED"] ? proc[p]["CACHED"] : 0)
                }
                printf "\n%-40s %8d %8d %8d\n", "TOTAL",
                    (status["COMPLETED"] ? status["COMPLETED"] : 0),
                    (status["FAILED"] ? status["FAILED"] : 0),
                    (status["CACHED"] ? status["CACHED"] : 0)
            }'
        echo ""

        echo "--- Failed processes (if any) ---"
        tail -n +2 "${TRACE}" | awk -F'\t' '$6 == "FAILED" {print $4, $5, $NF}' || echo "None"
        echo ""

        echo "--- Currently running ---"
        tail -n +2 "${TRACE}" | awk -F'\t' '$6 == "RUNNING" {print $4, $5}' || echo "None"
        echo ""
    else
        echo "Trace file not found: ${TRACE}"
        echo "Pipeline may not have started yet."
    fi

    # Check timeline
    if [[ -f "${OUTDIR}/pipeline_info/timeline.html" ]]; then
        echo "Timeline report: ${OUTDIR}/pipeline_info/timeline.html"
    fi
    if [[ -f "${OUTDIR}/pipeline_info/report.html" ]]; then
        echo "Execution report: ${OUTDIR}/pipeline_info/report.html"
    fi
    echo ""
}

show_disk() {
    echo "=== Disk Usage Summary ==="
    echo ""

    echo "--- Work directory ---"
    if [[ -d "${WORKDIR}" ]]; then
        du -sh "${WORKDIR}" 2>/dev/null || echo "Cannot read ${WORKDIR}"
        echo "  Subdirectories:"
        du -sh "${WORKDIR}"/*/ 2>/dev/null | sort -rh | head -10 || true
    else
        echo "Work directory not found: ${WORKDIR}"
    fi
    echo ""

    echo "--- Results directory ---"
    if [[ -d "${OUTDIR}" ]]; then
        du -sh "${OUTDIR}" 2>/dev/null || echo "Cannot read ${OUTDIR}"
        echo "  Subdirectories:"
        du -sh "${OUTDIR}"/*/ 2>/dev/null | sort -rh | head -10 || true
    else
        echo "Results directory not found: ${OUTDIR}"
    fi
    echo ""

    echo "--- FASTQ files ---"
    if [[ -d "${FASTQDIR}" ]]; then
        du -sh "${FASTQDIR}" 2>/dev/null
        ls -lh "${FASTQDIR}"/*.fastq.gz 2>/dev/null || true
    fi
    echo ""

    echo "--- References ---"
    if [[ -d "${REFDIR}" ]]; then
        du -sh "${REFDIR}" 2>/dev/null
    fi
    echo ""

    echo "--- Singularity cache ---"
    SINGCACHE="/rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/singularity_cache"
    if [[ -d "${SINGCACHE}" ]]; then
        du -sh "${SINGCACHE}" 2>/dev/null
        echo "  Container count: $(find "${SINGCACHE}" -name '*.sif' 2>/dev/null | wc -l)"
    fi
    echo ""

    echo "--- Ephemeral quota ---"
    df -h /rds/general/user/${HPC_USER:-YOUR_USERNAME}/ephemeral/ 2>/dev/null || \
        echo "Cannot check ephemeral quota"
    echo ""

    echo "--- Home quota ---"
    df -h /rds/general/user/${HPC_USER:-YOUR_USERNAME}/home/ 2>/dev/null || \
        quota -s 2>/dev/null || \
        echo "Cannot check home quota"
    echo ""
}

case "${CMD}" in
    queue|q)
        show_queue
        ;;
    log|l)
        show_log
        ;;
    progress|p)
        show_progress
        ;;
    disk|d)
        show_disk
        ;;
    all|a)
        show_queue
        show_progress
        show_disk
        show_log
        ;;
    watch|w)
        echo "Watching pipeline (refresh every 30s, Ctrl+C to stop)..."
        while true; do
            clear
            date
            echo ""
            show_queue
            show_progress
            sleep 30
        done
        ;;
    *)
        echo "Usage: $0 [queue|log|progress|disk|all|watch]"
        echo ""
        echo "  queue    (q)  - PBS job status"
        echo "  log      (l)  - Tail Nextflow log"
        echo "  progress (p)  - Pipeline process completion"
        echo "  disk     (d)  - Disk usage summary"
        echo "  all      (a)  - Show everything (default)"
        echo "  watch    (w)  - Auto-refresh queue + progress every 30s"
        ;;
esac

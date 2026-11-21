#!/bin/bash
BASE=/Volumes/RES/ws/paper1_ACLF_silence
IDX=$BASE/ref/salmon_index_GRCm39
OUT=$BASE/data/processed/salmon
mkdir -p "$OUT"
cd "$BASE/data/raw" || exit 1

for srr in SRR24030463 SRR24030461 SRR24030455 SRR24030459 SRR24030457 SRR24030453 SRR24030452 SRR24030462 SRR24030460 SRR24030454 SRR24030458 SRR24030456; do
  if [ -s "$OUT/$srr/quant.sf" ]; then
    echo "skip $srr (已完成)"
    continue
  fi
  if [ ! -s "${srr}_1_REAL.fastq.gz" ] || [ ! -s "${srr}_2_REAL.fastq.gz" ]; then
    echo "MISSING: $srr 的 fastq，跳过；补齐后重跑本脚本即可"
    continue
  fi
  echo "=== 定量中: $srr  $(date '+%H:%M:%S') ==="
  salmon quant -i "$IDX" -l A -p 8 --validateMappings \
    -1 "${srr}_1_REAL.fastq.gz" -2 "${srr}_2_REAL.fastq.gz" \
    -o "$OUT/$srr"
done
echo "=== 定量流程结束 ==="

import json
import pandas as pd
import re
import argparse

def extract_qc_with_trim(json_path, excel_output, sample_regex=None, depth_tsv=None):
    # 1. Load Data
    with open(json_path, 'r') as f:
        mqc_data = json.load(f)

    gen_stats = mqc_data.get('report_general_stats_data', [])
    raw_data = mqc_data.get('report_saved_raw_data', {})
    fastp_raw = raw_data.get('multiqc_fastp', {})
    fastqc_raw = raw_data.get('multiqc_fastqc', {})
    mosdepth_perchrom = raw_data.get('mosdepth_perchrom', {})

    merged = {}
    generic_metrics = {}

    # 2. Dynamic ID Extractor
    def get_id(s):
        if sample_regex:
            match = re.search(sample_regex, s)
            # Use group(1) if capturing group exists, otherwise use the whole match
            if match:
                return match.group(1) if match.groups() else match.group(0)
        return s  # If no regex provided or no match, assume the whole key is the ID

    # 3. Combine General Stats
    for entry in gen_stats:
        for s, metrics in entry.items():
            sid = get_id(s)

            # If a regex is provided, check if the key is a sample or generic
            is_sample = bool(re.search(sample_regex, s)) if sample_regex else True

            # Catch generic keys
            if sid == s and not is_sample:
                generic_metrics.update(metrics)
            else:
                if sid not in merged:
                    merged[sid] = {}
                merged[sid].update(metrics)

    # Merge generic metrics into our valid samples
    for sid in merged:
        merged[sid].update(generic_metrics)

    # 4. Extract before-filtering metrics from fastp raw data
    fastp_before = {}
    for sample_key, sample_data in fastp_raw.items():
        sid = get_id(sample_key)
        summary = sample_data.get('summary', {})
        before = summary.get('before_filtering', {})
        dup_rate = sample_data.get('duplication', {}).get('rate', 0)
        if before:
            fastp_before[sid] = {
                'before_gc_content': before.get('gc_content', 0),
                'before_q30_rate': before.get('q30_rate', 0),
                'before_read1_mean_length': before.get('read1_mean_length', 0),
                'fastp_dup_rate': dup_rate,
            }

    # 5. Extract FastQC %Duplicate before and after cleaning
    fastqc_dup = {}
    for fqc_key, fqc_data in fastqc_raw.items():
        sid = get_id(fqc_key)
        if sid not in fastqc_dup:
            fastqc_dup[sid] = {}

        pct_dup = 100.0 - fqc_data.get('total_deduplicated_percentage', 100.0)

        # Determine if this is before or after cleaning
        is_clean = fqc_key.endswith('_clean')
        is_r1 = '_R1' in fqc_key

        if is_clean:
            if is_r1:
                fastqc_dup[sid]['fastqc_dup_after_R1'] = pct_dup
            else:
                fastqc_dup[sid]['fastqc_dup_after_R2'] = pct_dup
        else:
            if is_r1:
                fastqc_dup[sid]['fastqc_dup_before_R1'] = pct_dup
            else:
                fastqc_dup[sid]['fastqc_dup_before_R2'] = pct_dup

    # 6. Load chr1-22 region depth (from external TSV or compute from mosdepth_perchrom)
    region_depth = {}
    if depth_tsv:
        try:
            df_depth = pd.read_csv(depth_tsv, sep='\t')
            for _, row in df_depth.iterrows():
                region_depth[row['Sample']] = {
                    'mean_region_depth': row['Mean_Region_Depth_chr1_22'],
                    'median_region_depth': row['Median_Region_Depth_chr1_22'],
                }
        except Exception as e:
            print(f"Warning: Could not read depth TSV '{depth_tsv}': {e}")

    if not region_depth and mosdepth_perchrom:
        import statistics
        for sample_key, chrom_data in mosdepth_perchrom.items():
            sid = get_id(sample_key)
            auto_vals = []
            for chrom, depth in chrom_data.items():
                if chrom.startswith('chr') and chrom[3:].isdigit():
                    chrom_num = int(chrom[3:])
                    if 1 <= chrom_num <= 22:
                        auto_vals.append(depth)
            if auto_vals:
                region_depth[sid] = {
                    'mean_region_depth': sum(auto_vals) / len(auto_vals),
                    'median_region_depth': statistics.median(auto_vals),
                }

    df_list = []
    for sid, dat in merged.items():
        # Filter out non-samples if a regex was provided
        if sample_regex and not re.search(sample_regex, sid):
            continue

        b_total = dat.get('before_filtering_total_reads', 0)
        a_total = dat.get('after_filtering_total_reads', 0)
        adp = dat.get('adapter_cutting_adapter_trimmed_reads', 0)

        # Get before-filtering metrics from fastp raw data
        bf = fastp_before.get(sid, {})
        before_gc = bf.get('before_gc_content', 0) * 100
        before_q30 = bf.get('before_q30_rate', 0) * 100
        before_length = bf.get('before_read1_mean_length', 0)
        fastp_dup = bf.get('fastp_dup_rate', 0) * 100

        # After-filtering metrics from general stats
        after_gc = dat.get('after_filtering_gc_content', 0) * 100
        after_q30 = dat.get('after_filtering_q30_rate', 0) * 100
        after_length = dat.get('after_filtering_read1_mean_length', 0)

        # FastQC %Duplicate (average of R1 and R2)
        fqc = fastqc_dup.get(sid, {})
        fqc_dup_before_r1 = fqc.get('fastqc_dup_before_R1', None)
        fqc_dup_before_r2 = fqc.get('fastqc_dup_before_R2', None)
        fqc_dup_after_r1 = fqc.get('fastqc_dup_after_R1', None)
        fqc_dup_after_r2 = fqc.get('fastqc_dup_after_R2', None)

        def avg_dup(r1, r2):
            vals = [v for v in [r1, r2] if v is not None]
            return sum(vals) / len(vals) if vals else 0

        fqc_dup_before = avg_dup(fqc_dup_before_r1, fqc_dup_before_r2)
        fqc_dup_after = avg_dup(fqc_dup_after_r1, fqc_dup_after_r2)

        # Median coverage from mosdepth
        median_cov = dat.get('median_coverage', 0)

        # Region depth (chr1-22)
        rd = region_depth.get(sid, {})

        row = {
            'Sample': sid,

            # BEFORE TRIM
            'Before_Total_Read(M)': b_total / 1e6 if b_total else 0,
            'Before_%GC': before_gc,
            'Before_Base_Quality(Q30)': before_q30,
            'Before_Length': before_length,
            'Before_%Dup(fastp)': fastp_dup,
            'Before_%Dup(FastQC-averageR1R2)': fqc_dup_before,

            # ADAPTER
            '%Adapter': (adp / b_total * 100) if b_total else 0,

            # AFTER TRIM
            'After_Total_Read(M)': a_total / 1e6 if a_total else 0,
            'After_%GC': after_gc,
            'After_Base_Quality(Q30)': after_q30,
            'After_Length': after_length,
            'After_%Dup(fastp)': dat.get('pct_duplication', 0),
            'After_%Dup(FastQC-averageR1R2)': fqc_dup_after,

            # ALIGNMENT & POST-TRIM METRICS
            'Map_Read (M)': dat.get('mapped_passed', 0) / 1e6,
            'Map%': dat.get('mapped_passed_pct', 0),
            '%Duplicate(samtools)': (dat.get('duplicates_passed', 0) / dat.get('flagstat_total', 1)) * 100 if 'flagstat_total' in dat else 0,
            'Mean_Cov(mosdepth)': dat.get('mean_coverage', 0),
            'Median_Cov(mosdepth)': median_cov,
            'Mean_Region_Depth_chr1-22': rd.get('mean_region_depth', 0),
            'Median_Region_Depth_chr1-22': rd.get('median_region_depth', 0),
            'Insert size': dat.get('summed_median', 0),
        }
        df_list.append(row)

    df = pd.DataFrame(df_list)

    if df.empty:
        print("Error: No data successfully matched and extracted. Check JSON keys or your regex pattern.")
        return

    # 7. Apply QC Filter
    def get_qc_status(row):
        m = row['Map%']
        c = row['Median_Cov(mosdepth)']
        i = row.get('Insert size', 0)
        d = row['%Duplicate(samtools)']
        bq = row['After_Base_Quality(Q30)']

        fails = []
        if pd.isna(m) or m < 90: fails.append("Low %Map")
        if pd.isna(c) or c < 10: fails.append("Low Cov")
        if pd.isna(i) or i < 100: fails.append("Small Ins")
        if pd.isna(d) or d > 20: fails.append("High Dup")
        if pd.isna(bq) or bq <= 90: fails.append("Low Base Qual")

        return "PASS" if not fails else f"FAIL ({', '.join(fails)})"

    df['Filter QC'] = df.apply(get_qc_status, axis=1)

    # 8. Save to Excel
    df.to_excel(excel_output, index=False)
    print(f"Generated complete report: {excel_output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract QC data from MultiQC JSON.")
    parser.add_argument("-i", "--input", required=True, help="Path to input MultiQC JSON file")
    parser.add_argument("-o", "--output", required=True, help="Path to output Excel file")
    parser.add_argument("-p", "--pattern", required=False, help=r"Regex pattern to extract sample names (e.g., '(.*?_Fx\w+M)')")
    parser.add_argument("-d", "--depth-tsv", required=False, help="Path to mosdepth region depth TSV (chr1-22 mean/median)")

    args = parser.parse_args()

    extract_qc_with_trim(args.input, args.output, args.pattern, args.depth_tsv)

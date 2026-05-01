import json
import pandas as pd
import re
import argparse
import statistics

def extract_qc_with_trim(json_path, excel_output, sample_regex=None):
    with open(json_path, 'r') as f:
        mqc_data = json.load(f)

    gen_stats = mqc_data.get('report_general_stats_data', [])
    raw_stats = mqc_data.get('report_saved_raw_data', {}).get('multiqc_general_stats', {})
    fastp_data = mqc_data.get('report_saved_raw_data', {}).get('multiqc_fastp', {})
    perchrom_data = mqc_data.get('report_saved_raw_data', {}).get('mosdepth_perchrom', {})

    def get_id(s):
        if sample_regex:
            match = re.search(sample_regex, s)
            if match:
                return match.group(1) if match.groups() else match.group(0)
        return s

    skip = re.compile(r'_(R1|R2)(_clean)?$|_flagstat$|insert_size|quality_yield', re.IGNORECASE)
    rows = {}

    # Pass 1: collect sample metrics
    for entry in gen_stats:
        for s, dat in entry.items():
            if skip.search(s) or (sample_regex and not re.search(sample_regex, s)):
                continue
            sid = get_id(s)
            b_total = dat.get('before_filtering_total_reads', 0)
            a_total = dat.get('after_filtering_total_reads', 0)
            adp = dat.get('adapter_cutting_adapter_trimmed_reads', 0)
            after_q30 = dat.get('after_filtering_q30_rate', 0) * 100
            rows[sid] = {
                'Sample': sid,
                'Before_Total_Read(M)': b_total / 1e6 if b_total else 0,
                'Before_%GC': dat.get('after_filtering_gc_content', 0) * 100,
                'Before_Base_Quality(Q30)': after_q30,
                'Before_Length': dat.get('after_filtering_read1_mean_length', 0),
                '%Adapter': (adp / b_total * 100) if b_total else 0,
                'After_Total_Read(M)': a_total / 1e6 if a_total else 0,
                'After_%GC': dat.get('after_filtering_gc_content', 0) * 100,
                'After_Base_Quality(Q30)': after_q30,
                'After_Length': dat.get('after_filtering_read1_mean_length', 0),
                'Map_Read (M)': 0,
                'Map%': 0,
                '%Duplicate': 0,
                'Base_Quality (Phred)': after_q30,
                'Median_Cov': 0,
                'Median_Cov_chr1_22': 0,
                'Insert size': 0,
            }

    # Pass 2: fill Map_Read and Map% from flagstat
    for entry in gen_stats:
        for s, dat in entry.items():
            if not s.endswith('_flagstat'):
                continue
            sid = get_id(s[:-len('_flagstat')])
            if sid not in rows:
                continue
            mapped = dat.get('mapped_passed', 0)
            a_total = rows[sid]['After_Total_Read(M)'] * 1e6
            rows[sid]['Map_Read (M)'] = mapped / 1e6
            rows[sid]['Map%'] = (mapped / a_total * 100) if a_total else 0

    autosome_re = re.compile(r'^chr([1-9]|1\d|2[0-2])$')

    # Pass 3: fill duplicate rate and coverage
    for s, data in raw_stats.items():
        sid = get_id(s)
        if sid in rows:
            rows[sid]['%Duplicate'] = data.get('fastp-pct_duplication', 0)
            rows[sid]['Median_Cov'] = data.get('mosdepth-median_coverage', 0)

    for s, data in perchrom_data.items():
        sid = get_id(s)
        if sid in rows:
            vals = [v for k, v in data.items() if autosome_re.match(k)]
            rows[sid]['Median_Cov_chr1_22'] = statistics.median(vals) if vals else 0

    # Pass 4: fill insert size and before length from fastp
    for s, data in fastp_data.items():
        sid = get_id(s)
        if sid in rows:
            rows[sid]['Insert size'] = data.get('insert_size', {}).get('peak', 0)
            rows[sid]['Before_Length'] = data.get('summary', {}).get('before_filtering', {}).get('read1_mean_length', 0)

    df_list = list(rows.values())
    #print(f"\n--- Collected {len(df_list)} rows ---")
    #for row in df_list:
    #    print(row)

    df = pd.DataFrame(df_list)
    if df.empty:
        print("Error: No data successfully matched and extracted. Check JSON keys or your regex pattern.")
        return

    def get_qc_status(row):
        fails = []
        if pd.isna(row['Map%']) or row['Map%'] < 90: fails.append("Low %Map")
        if pd.isna(row['Median_Cov']) or row['Median_Cov'] < 10: fails.append("Low Cov")
        if pd.isna(row['Insert size']) or row['Insert size'] < 100: fails.append("Small Ins")
        if pd.isna(row['%Duplicate']) or row['%Duplicate'] > 20: fails.append("High Dup")
        if pd.isna(row['Base_Quality (Phred)']) or row['Base_Quality (Phred)'] <= 90: fails.append("Low Base Qual")
        return "PASS" if not fails else f"FAIL ({', '.join(fails)})"

    df['Filter QC'] = df.apply(get_qc_status, axis=1)
    df.to_excel(excel_output, index=False)
    print(f"Generated complete report: {excel_output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract QC data from MultiQC JSON.")
    parser.add_argument("-i", "--input", required=True, help="Path to input MultiQC JSON file")
    parser.add_argument("-o", "--output", required=True, help="Path to output Excel file")
    parser.add_argument("-p", "--pattern", required=False, help=r"Regex pattern to extract sample names (e.g., '(.*?_Fx\w+M)')")
    args = parser.parse_args()
    extract_qc_with_trim(args.input, args.output, args.pattern)

import json
import pandas as pd
import re
import sys

def extract_qc_with_trim(json_path, excel_output):
    # 1. Load Data
    with open(json_path, 'r') as f:
        mqc_data = json.load(f)

    gen_stats = mqc_data.get('report_general_stats_data', [])
    raw_data = mqc_data.get('report_saved_raw_data', {})
    sources = mqc_data.get('report_data_sources', {})

    # 2. Build Name Map (To unify Fastp, Picard, Mosdepth, and Samtools)
    name_map = {}
    for tool, categories in sources.items():
        for cat, samples in categories.items():
            for internal_name, path in samples.items():
                # Extract base Sample ID (e.g., 20261125_WGS20X_OS00007_CF_Fx12M)
                match = re.search(r'(\d{8}.*?_Fx\w+M)', path)
                if match:
                    name_map[internal_name] = match.group(1)

    merged = {}
    def get_id(s):
        return name_map.get(s, s)

    # 3. Combine General Stats
    for entry in gen_stats:
        for s, metrics in entry.items():
            sid = get_id(s)
            if sid not in merged: merged[sid] = {}
            merged[sid].update(metrics)

    # 4. Extract Fastp Data (Before/After Trim)
    fastp = raw_data.get('multiqc_fastp', {})
    for s, metrics in fastp.items():
        sid = get_id(s)
        if sid not in merged: merged[sid] = {}
        merged[sid]['fastp_before'] = metrics.get('summary', {}).get('before_filtering', {})
        merged[sid]['fastp_after'] = metrics.get('summary', {}).get('after_filtering', {})
        merged[sid]['fastp_adapter'] = metrics.get('adapter_cutting', {}).get('adapter_trimmed_reads', 0)

    # 5. Extract Picard Quality Distribution (Average Phred)
    q_dist = raw_data.get('multiqc_picard_quality_score_distribution', {})
    for s, hist in q_dist.items():
        sid = get_id(s)
        if sid not in merged: merged[sid] = {}
        try:
            merged[sid]['phred'] = sum(float(v['QUALITY']) * float(v['COUNT_OF_Q']) for v in hist.values()) / sum(float(v['COUNT_OF_Q']) for v in hist.values())
        except:
            merged[sid]['phred'] = 0.0

    # 6. Build DataFrame
    df = pd.DataFrame()
    for sid, dat in merged.items():
        # Only process samples that have Fastp data
        if 'fastp_before' not in dat: continue 
        
        b = dat.get('fastp_before', {})
        a = dat.get('fastp_after', {})
        adp = dat.get('fastp_adapter', 0)
        
        row = {
            'Sample': sid,
            
            # BEFORE TRIM
            'Before_Total_Read(M)': b.get('total_reads', 0) / 1e6,
            'Before_%GC': b.get('gc_content', 0) * 100,
            'Before_Base_Quality(Q30)': b.get('q30_rate', 0) * 100,
            'Before_Length': b.get('read1_mean_length', 0),
            
            # ADAPTER
            '%Adapter': (adp / b.get('total_reads', 1)) * 100 if b.get('total_reads') else 0,
            
            # AFTER TRIM
            'After_Total_Read(M)': a.get('total_reads', 0) / 1e6,
            'After_%GC': a.get('gc_content', 0) * 100,
            'After_Base_Quality(Q30)': a.get('q30_rate', 0) * 100,
            'After_Length': a.get('read1_mean_length', 0),
            
            # ALIGNMENT & POST-TRIM METRICS
            'Map_Read (M)': dat.get('mapped_passed', 0) / 1e6,
            'Map%': dat.get('mapped_passed_pct', 0),
            '%Duplicate': (dat.get('duplicates_passed', 0) / dat.get('flagstat_total', 1)) * 100 if 'flagstat_total' in dat else 0,
            'Base_Quality (Phred)': dat.get('phred', 0),
            'Median_Cov': dat.get('median_coverage', 0),
            'Insert size': dat.get('summed_median', 0),
        }
        df = pd.concat([df, pd.DataFrame([row])])

    # 7. Apply QC Filter
    def get_qc_status(row):
        m, c, i, d = row['Map%'], row['Median_Cov'], row.get('Insert size', 0), row['%Duplicate']
        fails = []
        if m < 90: fails.append("Low %Map")
        if c < 10: fails.append("Low Cov")
        if i < 100: fails.append("Small Ins")
        if d > 20: fails.append("High Dup")
        return "PASS" if not fails else f"FAIL ({', '.join(fails)})"

    df['Filter QC'] = df.apply(get_qc_status, axis=1)
    
    # 8. Save to Excel
    df.to_excel(excel_output, index=False)
    print(f"Generated complete report: {excel_output}")


if __name__ == "__main__":
    # Input setting
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Execute the extraction
    extract_qc_with_trim(input_file, output_file)

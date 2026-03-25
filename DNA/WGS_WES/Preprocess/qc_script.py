import json
import pandas as pd
import re
import sys

def extract_qc_with_trim(json_path, excel_output):
    # 1. Load Data
    with open(json_path, 'r') as f:
        mqc_data = json.load(f)

    gen_stats = mqc_data.get('report_general_stats_data', [])
    
    merged = {}
    generic_metrics = {}
    
    # 2. Improved ID Extractor
    def get_id(s):
        match = re.search(r'(\d{8}.*?_Fx\w+M)', s)
        return match.group(1) if match else s

    # 3. Combine General Stats
    for entry in gen_stats:
        for s, metrics in entry.items():
            sid = get_id(s)
            
            # Catch generic keys like "insert_size" or "quality_yield"
            if sid == s and not re.search(r'Fx\w+M', s):
                generic_metrics.update(metrics)
            else:
                if sid not in merged: 
                    merged[sid] = {}
                merged[sid].update(metrics)
                
    # Merge generic metrics into our valid samples
    for sid in merged:
        merged[sid].update(generic_metrics)

    # 4. Build DataFrame
    df_list = []
    for sid, dat in merged.items():
        if not re.search(r'Fx\w+M', sid): 
            continue 
        
        b_total = dat.get('before_filtering_total_reads', 0)
        a_total = dat.get('after_filtering_total_reads', 0)
        adp = dat.get('adapter_cutting_adapter_trimmed_reads', 0)
        
        # Calculate Q30 once to use for both columns
        after_q30 = dat.get('after_filtering_q30_rate', 0) * 100
        
        row = {
            'Sample': sid,
            
            # BEFORE TRIM
            'Before_Total_Read(M)': b_total / 1e6 if b_total else 0,
            'Before_%GC': dat.get('after_filtering_gc_content', 0) * 100,
            'Before_Base_Quality(Q30)': dat.get('after_filtering_q30_rate', 0) * 100,
            'Before_Length': dat.get('after_filtering_read1_mean_length', 0),
            
            # ADAPTER
            '%Adapter': (adp / b_total * 100) if b_total else 0,
            
            # AFTER TRIM
            'After_Total_Read(M)': a_total / 1e6 if a_total else 0,
            'After_%GC': dat.get('after_filtering_gc_content', 0) * 100,
            'After_Base_Quality(Q30)': after_q30,
            'After_Length': dat.get('after_filtering_read1_mean_length', 0),
            
            # ALIGNMENT & POST-TRIM METRICS
            'Map_Read (M)': dat.get('mapped_passed', 0) / 1e6,
            'Map%': dat.get('mapped_passed_pct', 0),
            '%Duplicate': (dat.get('duplicates_passed', 0) / dat.get('flagstat_total', 1)) * 100 if 'flagstat_total' in dat else dat.get('pct_duplication', 0),
            
            # Copied from After Q30
            'Base_Quality (Phred)': after_q30, 
            
            'Median_Cov': dat.get('median_coverage', 0),
            'Insert size': dat.get('summed_median', 0),
        }
        df_list.append(row)

    df = pd.DataFrame(df_list)

    if df.empty:
        print("Error: No data successfully matched and extracted. Check JSON keys.")
        return

    # 5. Apply QC Filter
    def get_qc_status(row):
        m = row['Map%']
        c = row['Median_Cov']
        i = row.get('Insert size', 0)
        d = row['%Duplicate']
        bq = row['Base_Quality (Phred)']
        
        fails = []
        if pd.isna(m) or m < 90: fails.append("Low %Map")
        if pd.isna(c) or c < 10: fails.append("Low Cov")
        if pd.isna(i) or i < 100: fails.append("Small Ins")
        if pd.isna(d) or d > 20: fails.append("High Dup")
        if pd.isna(bq) or bq <= 90: fails.append("Low Base Qual") # Checks if Base Quality is > 90
        
        return "PASS" if not fails else f"FAIL ({', '.join(fails)})"

    df['Filter QC'] = df.apply(get_qc_status, axis=1)
    
    # 6. Save to Excel
    df.to_excel(excel_output, index=False)
    print(f"Generated complete report: {excel_output}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <input.json> <output.xlsx>")
        sys.exit(1)
        
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    extract_qc_with_trim(input_file, output_file)

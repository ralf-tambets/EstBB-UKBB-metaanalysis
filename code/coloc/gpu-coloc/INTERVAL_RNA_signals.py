import os
import pandas as pd
from tqdm import tqdm

signals_dir = "INTERVAL/signals"
summary_file_path = "INTERVAL/INTERVAL_signals.tsv"
os.makedirs(signals_dir, exist_ok=True)    

for ds in["INTERVAL_RNA_WGS_ge_blood", "INTERVAL_RNA_ge_blood"]:
    df = pd.read_parquet(f"interval/{ds}.lbf_variable.parquet")

    trait_data = df.groupby(df['molecular_trait_id'])
    
    for trait_id, group in tqdm(trait_data, desc=f"Processing molecular traits for {ds}"):
        for i in range(1, 11):
            lbf_column = f'lbf_variable{i}'
            if lbf_column in group.columns:
                df_filtered = group[['molecular_trait_id', 'region', 'variant', 'chromosome', 'position', lbf_column]].copy()
                df_filtered.rename(columns={lbf_column: 'lbf'}, inplace=True)
                
                if (df_filtered['lbf'] == 0).all():
                    continue
                
                signal_strength = df_filtered['lbf'].abs().max()

                if signal_strength < 5:
                    continue
                
                chromosome = df_filtered['chromosome'].iloc[0]
                location_min = df_filtered['position'].min()
                location_max = df_filtered['position'].max()

                signal = f"{ds}_{trait_id}_L{i}"
                output_file_name = f"{signal}.pickle"
                output_file_path = os.path.join(signals_dir, output_file_name)

                df_filtered['lbf'] = pd.to_numeric(df_filtered['lbf'], errors='coerce')

                mat1_df = pd.DataFrame(df_filtered.set_index('variant')['lbf']).T
                variant_id = mat1_df.T["lbf"].idxmax()

                mat1_df.to_pickle(output_file_path)

                summary_data = pd.DataFrame([{
                    'signal': signal,
                    'chromosome': chromosome,
                    'location_min': location_min,
                    'location_max': location_max,
                    'signal_strength': signal_strength,
                    'lead_variant': variant_id
                }])
                
                header_needed = not os.path.exists(summary_file_path)  
                summary_data.to_csv(summary_file_path, sep='\t', mode='a', header=header_needed, index=False)
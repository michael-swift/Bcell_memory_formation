import argparse
import pandas as pd
import scanpy as sc
from pathlib import Path
import logging

def setup_logging(log_level):
    logging.basicConfig(level=log_level,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

def shared_locations(group, subanatomical, loc1, loc2):
    if group.empty:
        return pd.NA
    if not subanatomical:
        return (loc1 in group['tissue'].values) and (loc2 in group['tissue'].values)
    else:
        return (loc1 in group['subanatomical_location'].values) and (loc2 in group['subanatomical_location'].values)

def main(args):
    setup_logging(args.log_level)
    logger = logging.getLogger(__name__)

    # Configuration
    data_dir = Path(args.data_dir)
    h5ad = data_dir / 'gex/bcells.h5ad.gz'
    vdj_infile = data_dir / 'vdj/integrated_cell_calls_ambient_annotated.tsv.gz'
    sharing_outfile = data_dir / 'annotation/sharing_labels_gex.tsv.gz'

    # Ensure output directory exists
    output_dir = sharing_outfile.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    seq_identifier = 'vdj_sequence'

    logger.info(f"Reading h5ad file: {h5ad}")
    adata = sc.read_h5ad(h5ad, backed='r')

    if args.restrict_to_TBd6:
        logger.info("Restricting analysis to TBd6 donor")
        adata = adata[adata.obs.donor == 'TBd6']
    adata = adata.to_memory()
    adata = adata[adata.obs.vdj_sequence != 'nan']
    adata.obs = adata.obs[['cb', 'vdj_sequence', 'celltypist', 'sample_uid']]

    if args.compute_sharing:
        logger.info(f"Computing sharing information from {vdj_infile}")
        df = pd.read_table(vdj_infile, index_col=0, usecols=["vdj_sequence", "lineage_id", 'v_mismatch', 'donor', 'probable_hq_single_b_cell', 'subanatomical_location', 'tissue', 'locus', 'c_call', 'sample_uid', 'cb'])
        df = df.dropna(subset=['locus'])
        df_LN = df[df.tissue == 'LN']
        logger.info("Subanatomical location value counts:")
        logger.info(df_LN.subanatomical_location.value_counts())

        shared_LNs = df_LN.groupby(seq_identifier)['subanatomical_location'].nunique() > 1
        df['shared_LN_LN'] = df[seq_identifier].map(shared_LNs)

        found_single_LN = df_LN.groupby(seq_identifier)['subanatomical_location'].nunique() == 1
        df['found_in_single_LN'] = df[seq_identifier].map(found_single_LN)

        shared_categories = {
            'shared_LN_SP': {'subanatomical': False, 'loc1': 'LN', 'loc2': 'SP'},
            'shared_SP_PB': {'subanatomical': False, 'loc1': 'SP', 'loc2': 'PB'},
            'shared_BM_PB': {'subanatomical': False, 'loc1': 'BM', 'loc2': 'PB'},
            'shared_LN_PB': {'subanatomical': False, 'loc1': 'LN', 'loc2': 'PB'},
            'shared_LN_BM': {'subanatomical': False, 'loc1': 'LN', 'loc2': 'BM'},
            'shared_SP_BM': {'subanatomical': False, 'loc1': 'SP', 'loc2': 'BM'},
        }
        for category, params in shared_categories.items():
            logger.info(f"Computing {category}")
            shared = df.groupby(seq_identifier).apply(shared_locations, params['subanatomical'], params['loc1'], params['loc2'])
            df[category] = df[seq_identifier].map(shared)

        logger.info(f"Saving sharing information to {sharing_outfile}")
        df.reset_index().to_csv(sharing_outfile, sep="\t", index=False)
    
    logger.info(f"Reading sharing information from {sharing_outfile}")
    shared = pd.read_table(sharing_outfile, index_col=0)

    logger.info(f"{shared.shape[0]} antibody assemblies to analyze for sharing")
    shared = shared[shared.probable_hq_single_b_cell == True]
    logger.info(f"{shared.shape[0]} number of hq single B cells transcriptomes to detect GEX differences amongst shared cells")

    logger.info("Merging data")
    merged_df = pd.merge(adata.obs, shared, left_on=['cb', 'vdj_sequence'], right_on=['cb', 'vdj_sequence'], suffixes=('_left', ''), how='left')
    cols_to_remove = [col for col in merged_df.columns if '_left' in col]
    merged_df.drop(columns=cols_to_remove, inplace=True)
    merged_df.index = adata.obs.index

    logger.info("Converting column types")
    for col in merged_df.columns:
        if merged_df[col].dtype == 'object':
            merged_df[col] = merged_df[col].astype('category')
        elif pd.api.types.is_float_dtype(merged_df[col]):
            merged_df[col] = merged_df[col].astype('float32')
        elif pd.api.types.is_integer_dtype(merged_df[col]):
            merged_df[col] = merged_df[col].astype('int32')
        else:
            merged_df[col] = merged_df[col].astype('category')

    adata = adata.copy()
    adata.obs = merged_df
    adata = adata[~adata.obs.sample_uid.isna()]

    output_dir = data_dir / 'annotation'
    output_dir.mkdir(parents=True, exist_ok=True)

    output_filename = "TBd6_sharing.h5ad.gz" if args.restrict_to_TBd6 else "all_sharing.h5ad.gz"
    output_file = output_dir / "sharing" / output_filename

    logger.info(f"Writing output to {output_file}")
    adata.write_h5ad(output_file, compression='gzip')
    logger.info(f"Output written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process B cell data and compute sharing information.")
    parser.add_argument("--data_dir", type=str, default="../../data/", help="Path to the data directory")
    parser.add_argument("--restrict_to_TBd6", action="store_true", help="Restrict analysis to TBd6 donor")
    parser.add_argument("--compute_sharing", action="store_true", help="Compute sharing information")
    parser.add_argument("--log_level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="Set the logging level")
    args = parser.parse_args()

    main(args)
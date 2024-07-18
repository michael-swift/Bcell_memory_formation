import pathlib
import pandas as pd
import argparse

def generate_latex_table(output_path, sample_relationships):
    relationships = pd.read_table(sample_relationships)
    relationships.loc[:, 'negative_selection'] = (~relationships['donor_tissue'].str.contains("LN") & (~relationships['sample_uid'].str.contains("TBd5_frozen_PB")))
    relationships['frozen'] = relationships.cell_suspension.str.contains("frozen")
    relationships = relationships.sort_values(['donor_tissue', 'frozen'])
    relationships = relationships[['cell_suspension', 'sample_uid', 'encapsulation_sibling_sample', 'is_gex_sample', 'negative_selection']]
    relationships.columns = ["Cell Suspension", "Unique Sample ID", "Technical Sibling", "GEX", "B cell enriched"]
    relationships['Technical Sibling'] = relationships['Technical Sibling'].fillna("")

    df = relationships

    # Initialize the LaTeX table string
    latex_table = "\\floatsetup[longtable]{LTcapwidth=table}\n"
    latex_table += "\\begin{longtable}{lllp{1cm}p{1cm}p{1.3cm}}\n"
    latex_table += "\\midrule\n"
    latex_table += "Cell Suspension & Unique Sample ID & Technical Sibling & GEX & B{\\color{white}B}cell enriched \\\\\n"
    latex_table += "\\midrule\n"
    latex_table += "\\endhead\n"

    # Loop through each row to populate the LaTeX table
    for index, row in df.iterrows():
        gex_check = "\\checkmark" if row['GEX'] else ""
        b_cell_check = "\\checkmark" if row['B cell enriched'] else ""
        latex_table += f"\\verb|{row['Cell Suspension']}| & \\verb|{row['Unique Sample ID']}| & \\verb|{row['Technical Sibling']}| & {gex_check} & {b_cell_check} \\\\\n"

    # Write the table to a text file
    with open(output_path, 'w') as file:
        file.write(latex_table)

def main():
    parser = argparse.ArgumentParser(description="Generate LaTeX table from sample relationships data.")
    parser.add_argument("--sample_relationships", type=str, help="Path to sample relationships data.", default="../../snakemake_workflow/samplesheets/sample_relationships.tsv")
    parser.add_argument("--output_path", type=str, help="Path to output the LaTeX table text file.", default="sample_information.txt")
    args = parser.parse_args()
    generate_latex_table(args.output_path, args.sample_relationships)

if __name__ == "__main__":
    main()
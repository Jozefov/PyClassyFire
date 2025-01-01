# import click
# import json
# import os
# from .src.batch import Job, Scheduler, process_batches_with_saving_and_retry
# from .src.utils import MoleCule, chunk_tasks, load_existing_results, save_intermediate_results
# from requests.exceptions import HTTPError
#
#
# @click.command()
# @click.option('--identifier', type=click.Choice(['SMILES', 'InChI']), default='SMILES', help='Type of identifier')
# @click.argument('input_file', type=click.Path(exists=True))
# @click.argument('output_file', type=click.Path())
# def main(identifier, input_file, output_file):
#     """
#     Classify molecules using ClassyFire API with resumption and enhanced error handling.
#
#     INPUT_FILE: Path to the input file (JSON or TXT/CSV) containing SMILES or InChI.
#     OUTPUT_FILE: Path to the output file (JSON or TSV) to save classification results.
#     """
#     # Determine input file type
#     input_suffix = os.path.splitext(input_file)[1].lower()
#     with open(input_file, 'r') as f:
#         if input_suffix == '.json':
#             identifiers = json.load(f)
#         elif input_suffix in ['.txt', '.csv']:
#             identifiers = [line.strip() for line in f]
#         else:
#             raise ValueError('Input file must be a JSON or TXT/CSV file.')
#
#     # Canonicalize identifiers
#     if identifier == 'SMILES':
#         canonical_ids = [MoleCule.from_smiles(x).canonical_smiles for x in identifiers if x]
#     else:
#         canonical_ids = [MoleCule.from_inchi(x).canonical_smiles for x in identifiers if x]
#
#     # Determine batch size based on identifier type
#     current_batch_size = 100 if identifier == 'SMILES' else 25
#
#     # Define output directory for intermediate results
#     intermediate_dir = os.path.join(os.path.dirname(output_file), 'intermediate_results')
#
#     # Process batches with resumption and error handling
#     intermediate_files = process_batches_with_saving_and_retry(
#         smiles_list=canonical_ids,
#         batch_size=current_batch_size,
#         save_interval=20,  # Adjust as needed
#         output_dir=intermediate_dir,
#         max_retries=3,
#         retry_delay=300
#     )
#
#     # Merge intermediate results
#     merged_results = {}
#     for file in intermediate_files:
#         try:
#             with open(file, 'r') as f:
#                 data = json.load(f)
#                 merged_results.update(data)
#             click.echo(f"Successfully merged results from {file}")
#         except Exception as e:
#             click.echo(f"Error merging results from {file}: {e}")
#
#     # Map back to original identifiers
#     smiles2identifier = dict(zip([MoleCule.from_smiles(
#         x).canonical_smiles if identifier == 'SMILES' else MoleCule.from_inchi(x).canonical_smiles for x in identifiers
#                                   if x],
#                                  identifiers))
#     mapped_results = {smiles2identifier.get(k, k): v for k, v in merged_results.items()}
#
#     # Export results
#     output_suffix = os.path.splitext(output_file)[1].lower()
#     if output_suffix == '.json':
#         with open(output_file, 'w') as f_out:
#             json.dump(mapped_results, f_out, indent=4)
#     elif output_suffix == '.tsv':
#         with open(output_file, 'w') as f_out:
#             f_out.write('Identifier\tSuperclass\tClass\tSubclass\n')
#             for k, v in mapped_results.items():
#                 superclass = v.get('superclass', 'Unknown')
#                 cls = v.get('class', 'Unknown')
#                 subclass = v.get('subclass', 'Unknown')
#                 f_out.write(f"{k}\t{superclass}\t{cls}\t{subclass}\n")
#     else:
#         raise ValueError('Output file must be a JSON or TSV file.')
#
#     click.echo(f"Annotated data has been saved to {output_file}")
#
#
# if __name__ == '__main__':
#     main()
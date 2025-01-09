import click
import json
import pandas as pd
import os

from .src.utils import (
    MoleCule,
    merge_intermediate_files,
    check_all_smiles_present,
    save_smiles_mapping,

)
from .src.batch import process_batches_with_saving_and_retry

@click.group()
def cli():
    """PyClassyFire: A tool to classify SMILES using the ClassyFire API."""
    pass

@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path())
@click.option('--batch_size', default=100, show_default=True, help='Number of SMILES per batch (max 100).')
@click.option('--max_retries', default=3, show_default=True, help='Maximum number of retries for failed batches.')
@click.option('--retry_delay', default=10, show_default=True, help='Delay between retries in seconds.')
def classify(input_file, output_dir, batch_size, max_retries, retry_delay):
    """
    Classify SMILES using the ClassyFire API.

    INPUT_FILE: Path to the input file containing SMILES strings (formats: txt, tsv, csv, json).
    OUTPUT_DIR: Path to the output directory where results and logs will be saved.
    """
    # 1. Determine Input File Type and Load SMILES
    input_suffix = os.path.splitext(input_file)[1].lower()
    if input_suffix == '.json':
        with open(input_file, 'r') as f:
            identifiers = json.load(f)
            # Handle different JSON structures
            if isinstance(identifiers, dict):
                # If JSON is a dictionary, extract all values
                identifiers = list(identifiers.values())
            elif not isinstance(identifiers, list):
                raise ValueError("JSON input must be a list or a dictionary of SMILES strings.")
    elif input_suffix in ['.txt', '.tsv', '.csv']:
        if input_suffix == '.csv':
            df = pd.read_csv(input_file, header=None)
        else:  # .txt or .tsv
            df = pd.read_csv(input_file, sep='\t', header=None)
        identifiers = df.iloc[:, 0].tolist()
    else:
        raise ValueError("Input file must be a JSON, TXT, TSV, or CSV file containing SMILES strings.")

    # 2. Canonicalize SMILES and Remove Invalid Entries
    smiles_list = []
    invalid_smiles = 0
    for smi in identifiers:
        smi = smi.strip()
        if not smi:
            continue  # Skip empty lines
        molecule = MoleCule.from_smiles(smi)
        if molecule.mol:
            canonical_smi = molecule.canonical_smiles
            smiles_list.append(canonical_smi)
        else:
            invalid_smiles += 1
            click.echo(f"Invalid SMILES skipped: {smi}")

    # 3. Remove duplicates
    smiles_list = list(set(smiles_list))
    click.echo(f'Total unique valid SMILES to process: {len(smiles_list)}')

    if invalid_smiles > 0:
        click.echo(f"Number of invalid SMILES skipped: {invalid_smiles}")

    # 4. Define Paths for Intermediate Results and Final Output
    intermediate_dir = os.path.join(output_dir, 'intermediate_results')
    os.makedirs(intermediate_dir, exist_ok=True)

    final_output_path = os.path.join(output_dir, 'output.json')

    # 5. Save the Mapping from Original to Canonical SMILES
    # Deduplicate original_smiles_list using set to avoid redundancy
    original_smiles_list = list(set([smi.strip() for smi in identifiers if smi.strip()]))
    save_smiles_mapping(original_smiles_list, output_dir)

    # 6. Process Batches with Resumption and Retry Logic
    intermediate_files = process_batches_with_saving_and_retry(
        smiles_list=smiles_list,
        batch_size=batch_size,
        output_dir=intermediate_dir,
        max_retries=max_retries,
        retry_delay=retry_delay
    )

    # 7. Merge Intermediate Files into Final Output JSON
    merge_intermediate_files(intermediate_dir, final_output_path)

    # 8. Check Completeness of Processing
    # Load the mapping.json to map original SMILES to canonical SMILES
    mapping_file_path = os.path.join(output_dir, 'mapping.json')
    try:
        with open(mapping_file_path, 'r') as f:
            mapping = json.load(f)
    except Exception as e:
        raise ValueError(f"Failed to load mapping.json: {e}")

    # Prepare the list of canonical SMILES that correspond to the original SMILES
    # Only include SMILES that were successfully mapped
    mapped_canonical_smiles = [mapping[smi] for smi in original_smiles_list if
                               smi in mapping and mapping[smi] is not None]

    # Remove duplicates to match processing
    mapped_canonical_smiles = list(set(mapped_canonical_smiles))

    # Check if all mapped canonical SMILES are present in the final output
    missing_smiles_mapping = check_all_smiles_present(final_output_path, mapped_canonical_smiles)

    # 7. Save Missing SMILES Mapping if Any
    if missing_smiles_mapping:
        missing_smiles_path = os.path.join(output_dir, 'missing_smiles.json')
        try:
            with open(missing_smiles_path, 'w') as f:
                json.dump(missing_smiles_mapping, f, indent=4)
            click.echo(f"Missing SMILES mapping has been saved to {missing_smiles_path}.")
            click.echo("The 'missing_smiles.json' contains a mapping from original SMILES to their canonical forms that were not processed.")
        except Exception as e:
            click.echo(f"Error saving missing SMILES mapping: {e}")
    else:
        click.echo("No missing SMILES. All SMILES have been successfully processed.")

    click.echo(f"Final classification results have been saved to {final_output_path}.")


@cli.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_path', type=click.Path())
def map(input_file, output_path):
    """
    Generate a mapping from original SMILES to canonical SMILES without sending to the API.

    INPUT_FILE: Path to the input file containing SMILES strings (formats: txt, tsv, csv, json).
    OUTPUT_PATH: Path to the output directory or JSON file where mapping.json will be saved.
    """
    # 1. Determine Input File Type and Load SMILES
    input_suffix = os.path.splitext(input_file)[1].lower()
    if input_suffix == '.json':
        with open(input_file, 'r') as f:
            identifiers = json.load(f)
            # Handle different JSON structures
            if isinstance(identifiers, dict):
                # If JSON is a dictionary, extract all values
                identifiers = list(identifiers.values())
            elif not isinstance(identifiers, list):
                raise ValueError("JSON input must be a list or a dictionary of SMILES strings.")
    elif input_suffix in ['.txt', '.tsv', '.csv']:
        if input_suffix == '.csv':
            df = pd.read_csv(input_file, header=None)
        else:  # .txt or .tsv
            df = pd.read_csv(input_file, sep='\t', header=None)
        identifiers = df.iloc[:, 0].tolist()
    else:
        raise ValueError("Input file must be a JSON, TXT, TSV, or CSV file containing SMILES strings.")

    # 2. Prepare Original SMILES List (deduplicated)
    original_smiles_list = list(set([smi.strip() for smi in identifiers if smi.strip()]))
    click.echo(f'Total unique original SMILES: {len(original_smiles_list)}')

    # 3. Save the Mapping from Original to Canonical SMILES
    save_smiles_mapping(original_smiles_list, output_path)


@cli.command()
@click.argument('intermediate_dir', type=click.Path(exists=True, file_okay=False))
@click.argument('final_output_path', type=click.Path())
def merge(intermediate_dir, final_output_path):
    """
    Merge all intermediate JSON files into a single final JSON file.

    INTERMEDIATE_DIR: Path to the directory containing intermediate JSON files.
    FINAL_OUTPUT_PATH: Path to the output JSON file where merged results will be saved.
    """
    # 1. Merge Intermediate Files
    merge_intermediate_files(intermediate_dir, final_output_path)


if __name__ == '__main__':
    cli()
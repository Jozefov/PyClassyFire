from rdkit import Chem
import os
import json
import logging
import re
from typing import List, Dict
from rdkit import RDLogger

# Suppress RDKit warnings globally
RDLogger.DisableLog('rdApp.*')

class MoleCule:
    def __init__(self, mol):
        self.mol = mol
        # self.ts = rdMolStandardize.TautomerEnumerator()

    @classmethod
    def from_smiles(cls, smiles: str):
        mol = Chem.MolFromSmiles(smiles)
        return cls(mol)
    
    @classmethod
    def from_inchi(cls, inchi: str):
        mol = Chem.MolFromInchi(inchi)
        return cls(mol)

    @property
    def canonical_smiles(self):
        # mol = self.ts.Canonicalize(self.mol)
        # # Canonicalize SMILES with stereochemistry
        # Chem.SanitizeMol(mol, Chem.SANITIZE_SETAROMATICITY)
        return Chem.MolToSmiles(self.mol, isomericSmiles=False, canonical=True)

    def __str__(self):
        return self.canonical_smiles
    

def chunk_tasks(data: list, size: int) -> list[list]:
    """
    Split a list into chunks of a given size
    :param data: The list to split
    :param size: The size of each chunk
    :return: A list of chunks
    """
    return [data[i:i+size] for i in range(0, len(data), size)]
    


def take_class(raw: dict | None) -> str:
    """
    Extract the ClassyFire class from the raw JSON response
    :param raw: The raw JSON response
    :return: The ClassyFire class
    """
    try:
        return raw['name']
    except (KeyError, TypeError):
        return 'Unknown'


def load_existing_results(output_dir: str) -> Tuple[set, int]:
    """
    Loads existing intermediate JSON files and returns a set of already processed original SMILES.
    Also returns the highest batch number to continue numbering correctly.

    Parameters:
    - output_dir (str): Directory where intermediate JSON files are stored.

    Returns:
    - tuple:
        - already_processed (set): Set of original SMILES strings that have been processed.
        - max_batch_num (int): The highest batch number found in existing files.
    """
    already_processed = set()
    max_batch_num = 0
    if not os.path.exists(output_dir):
        return already_processed, max_batch_num

    batch_id_pattern = re.compile(r'intermediate_(\d+)\.json$')

    for file in os.listdir(output_dir):
        match = batch_id_pattern.match(file)
        if match:
            batch_num = int(match.group(1))
            if batch_num > max_batch_num:
                max_batch_num = batch_num
            file_path = os.path.join(output_dir, file)
            try:
                with open(file_path, 'r') as f:
                    data = json.load(f)
                    # Each file has a single batch_id key mapping to a list of molecule dicts
                    for batch_id, molecules in data.items():
                        for molecule in molecules:
                            original_smiles = molecule.get('original_smiles')
                            if original_smiles:
                                already_processed.add(original_smiles)
                logging.info(f"Loaded results from {file_path}")
            except Exception as e:
                logging.error(f"Failed to load {file_path}: {e}")
                print(f"Failed to load {file_path}: {e}")
    return already_processed, max_batch_num


def save_intermediate_results(batch_num, molecules, original_smiles, output_dir):
    """
    Saves the results of a single batch to an intermediate JSON file.

    Each molecule dictionary will include an 'original_smiles' key corresponding to the SMILES sent for classification.

    Parameters:
    - batch_num (int): Unique sequential identifier for the batch.
    - molecules (list): List of molecule dictionaries with classification details.
    - original_smiles (list): List of SMILES strings sent for classification.
    - output_dir (str): Directory to save intermediate JSON files.
    """
    file_name = f'intermediate_{batch_num}.json'
    file_path = os.path.join(output_dir, file_name)

    augmented_molecules = []
    for mol, orig_smi in zip(molecules, original_smiles):
        mol['original_smiles'] = orig_smi
        # Reorder dictionary to have 'original_smiles' first
        mol = {'original_smiles': mol['original_smiles'], **{k: v for k, v in mol.items() if k != 'original_smiles'}}
        augmented_molecules.append(mol)

    try:
        with open(file_path, 'w') as f:
            json.dump({f'Batch_{batch_num}': augmented_molecules}, f, indent=4)
        logging.info(f"Saved intermediate results to {file_path}")
        print(f"Saved intermediate results to {file_path}")
    except Exception as e:
        logging.error(f"Failed to save {file_path}: {e}")
        print(f"Failed to save {file_path}: {e}")


# def extract_smiles_classification(molecules):
#     """
#     Converts list of molecule dicts to a SMILES: classification_details dict.
#
#     Parameters:
#     - molecules (list): List of molecule dictionaries.
#
#     Returns:
#     - dict: Mapping from SMILES to their classification details.
#     """
#     extracted = {}
#     for molecule in molecules:
#         smiles = molecule.get('smiles')
#         if smiles:
#             classification = {
#                 "superclass": molecule.get('superclass', {}).get('name', 'Unknown'),
#                 "class": molecule.get('class', {}).get('name', 'Unknown'),
#                 "subclass": molecule.get('subclass', {}).get('name', 'Unknown')
#             }
#             extracted[smiles] = classification
#     return extracted


def merge_intermediate_files(
    output_dir: str,
    final_output_path: str
) -> None:
    """
    Merges all intermediate JSON files in the specified directory into a single JSON file.
    Preserves the 'original_smiles' field for each molecule.

    Parameters:
    - output_dir (str): Directory containing intermediate JSON files.
    - final_output_path (str): Path to save the merged final JSON file.

    Returns:
    - None
    """
    # Compile regex pattern to match intermediate files
    pattern = re.compile(r'intermediate_(\d+)\.json$')

    # List to hold all molecule objects
    all_molecules = []

    # Find and sort all matching intermediate files
    intermediate_files = []
    for file in os.listdir(output_dir):
        match = pattern.match(file)
        if match:
            batch_num = int(match.group(1))
            intermediate_files.append((batch_num, os.path.join(output_dir, file)))

    # Sort files based on batch number in ascending order
    intermediate_files.sort(key=lambda x: x[0])

    if not intermediate_files:
        print("No intermediate files found to merge.")
        return

    # Iterate through each intermediate file and collect molecule objects
    for batch_num, file_path in intermediate_files:
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
                # Extract the list of molecules
                for key, molecules in data.items():
                    if isinstance(molecules, list):
                        all_molecules.extend(molecules)
                    else:
                        print(f"Unexpected format in {file_path}: Expected a list of molecules.")
        except json.JSONDecodeError as e:
            print(f"JSON decode error in {file_path}: {e}")
        except Exception as e:
            print(f"Error reading {file_path}: {e}")

    # Save the merged list of molecules to the final output JSON file
    try:
        with open(final_output_path, 'w') as f_out:
            json.dump(all_molecules, f_out, indent=4)
        print(f"Successfully merged {len(intermediate_files)} files into {final_output_path}.")
    except Exception as e:
        print(f"Failed to save merged JSON file {final_output_path}: {e}")


def check_all_smiles_present(
        final_output_path: str,
        original_smiles_list: List[str]
) -> Dict[str, str]:
    """
    Checks whether all SMILES from the original list are present in the final output JSON file.
    Returns a mapping from original SMILES to canonical SMILES for the missing SMILES.

    Parameters:
    - final_output_path (str): Path to the merged final JSON file.
    - original_smiles_list (List[str]): Original list of SMILES strings that were classified.

    Returns:
    - Dict[str, str]: Mapping from original SMILES to canonical SMILES for missing SMILES.
    """
    # Load the final output JSON
    try:
        with open(final_output_path, 'r') as f:
            molecules = json.load(f)
    except FileNotFoundError:
        print(f"Final output file {final_output_path} not found.")
        return {}
    except json.JSONDecodeError as e:
        print(f"JSON decode error in {final_output_path}: {e}")
        return {}
    except Exception as e:
        print(f"Error reading {final_output_path}: {e}")
        return {}

    # Extract original SMILES from the final output
    output_original_smiles = set()
    for molecule in molecules:
        original_smiles = molecule.get('original_smiles')
        if original_smiles:
            output_original_smiles.add(original_smiles)

    # Create set from original_smiles_list
    original_smiles_set = set([smi.strip() for smi in original_smiles_list if smi.strip()])

    # Identify missing original SMILES
    missing_original_smiles = original_smiles_set - output_original_smiles

    # Map original SMILES to their canonical SMILES
    missing_smiles_mapping = {}
    for smi in missing_original_smiles:
        canonical_smi = MoleCule.from_smiles(smi).canonical_smiles
        if canonical_smi:
            missing_smiles_mapping[smi] = canonical_smi

    if not missing_smiles_mapping:
        print("All SMILES are present in the final output.")
    else:
        print(f"Missing {len(missing_smiles_mapping)} SMILES in the final output:")
        print("Can be caused by server ERROR, try to lower batch and rerun.")
        print("It keeps already processed molecules and will try to retrieve only missing smiles.")

    return missing_smiles_mapping

from rdkit import Chem
import os
import json
import logging
import re

class MoleCule:
    def __init__(self, mol):
        self.mol = mol

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
        return Chem.MolToSmiles(self.mol, canonical=True)
    
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


def load_existing_results(output_dir):
    """
    Loads existing intermediate JSON files and returns a set of already processed SMILES.
    Also returns the highest batch number to continue numbering correctly.

    Parameters:
    - output_dir (str): Directory where intermediate JSON files are stored.

    Returns:
    - tuple:
        - already_processed (set): Set of SMILES strings that have been processed.
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
                            smiles = molecule.get('smiles')
                            if smiles:
                                already_processed.add(smiles)
                logging.info(f"Loaded results from {file_path}")
            except Exception as e:
                logging.error(f"Failed to load {file_path}: {e}")
                print(f"Failed to load {file_path}: {e}")
    return already_processed, max_batch_num


def save_intermediate_results(batch_num, molecules, output_dir):
    """
    Saves the results of a single batch to an intermediate JSON file.

    Parameters:
    - batch_num (int): Unique sequential identifier for the batch.
    - molecules (list): List of molecule dictionaries with classification details.
    - output_dir (str): Directory to save intermediate JSON files.
    """
    file_name = f'intermediate_{batch_num}.json'
    file_path = os.path.join(output_dir, file_name)
    try:
        with open(file_path, 'w') as f:
            json.dump({f'Batch_{batch_num}': molecules}, f, indent=4)
        logging.info(f"Saved intermediate results to {file_path}")
        print(f"Saved intermediate results to {file_path}")
    except Exception as e:
        logging.error(f"Failed to save {file_path}: {e}")
        print(f"Failed to save {file_path}: {e}")


def extract_smiles_classification(molecules):
    """
    Converts list of molecule dicts to a SMILES: classification_details dict.

    Parameters:
    - molecules (list): List of molecule dictionaries.

    Returns:
    - dict: Mapping from SMILES to their classification details.
    """
    extracted = {}
    for molecule in molecules:
        smiles = molecule.get('smiles')
        if smiles:
            classification = {
                "superclass": molecule.get('superclass', {}).get('name', 'Unknown'),
                "class": molecule.get('class', {}).get('name', 'Unknown'),
                "subclass": molecule.get('subclass', {}).get('name', 'Unknown')
            }
            extracted[smiles] = classification
    return extracted
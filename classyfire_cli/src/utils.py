from rdkit import Chem
import os
import json
import logging

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


def load_existing_results(output_dir):
    """
    Loads existing intermediate JSON files and returns a dictionary of already processed identifiers.
    """
    merged_results = {}
    if not os.path.exists(output_dir):
        return merged_results
    for file in os.listdir(output_dir):
        if file.startswith('intermediate_') and file.endswith('.json'):
            file_path = os.path.join(output_dir, file)
            try:
                with open(file_path, 'r') as f:
                    data = json.load(f)
                    merged_results.update(data)
                logging.info(f"Loaded results from {file_path}")
            except Exception as e:
                logging.error(f"Failed to load {file_path}: {e}")
    return merged_results

def save_intermediate_results(results, output_dir, batch_num=None):
    """
    Saves the merged results to an intermediate JSON file.
    """
    if batch_num:
        filename = f'intermediate_{batch_num}.json'
    else:
        filename = 'intermediate_final.json'
    file_path = os.path.join(output_dir, filename)
    try:
        with open(file_path, 'w') as f:
            json.dump(results, f, indent=4)
        logging.info(f"Saved intermediate results to {file_path}")
    except Exception as e:
        logging.error(f"Failed to save {file_path}: {e}")
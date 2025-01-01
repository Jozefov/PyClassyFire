from rdkit import Chem

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
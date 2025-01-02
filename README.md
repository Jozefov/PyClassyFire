# PyClassyFire

A Python client for the [ClassyFire API](http://classyfire.wishartlab.com/) for large-scale chemical compound classification.

## Introduction

PyClassyFire is a Python client designed to interact with the [ClassyFire API](http://classyfire.wishartlab.com/) for the large-scale classification of chemical compounds. It streamlines the process of submitting SMILES strings to the API, handling canonicalization, batch processing, and result retrieval.

## Features

- **Batch Processing:** Process large lists of SMILES strings.
- **Canonicalization:** Ensures consistent SMILES representation using RDKit.
- **Error Handling:** Handling of parsing and API errors.
- **Resuming Processed molecules** Enable to start from where we stopped


## Installation

### Using Conda

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/Jozefov/PyClassyFire.git
   cd PyClassyFire
	```
2.  **Create and Activate the Conda Environment:**
	```bash
	conda env create -f environment.yml
	conda activate classyfire_env
	```
 3. **Install the Package:**
	```bash
	pip install .
	```
4. **Verify Installation:**
    ```bash
	pyclassyfire --help
	```
    You should see the help message detailing usage options.

**Usage**
PyClassyFire provides a command-line interface (CLI) to interact with the ClassyFire API.

**Running the Script**

After installation, you can use the classyfire command to classify your chemical compounds.

```bash
	pyclassyfire <input_file> <output_dir> [OPTIONS]
```

**Parameters**

•  <input_file>: Path to the input file containing SMILES strings (formats: txt, tsv, csv, json).

•  <output_dir>: Path to the output directory where results and logs will be saved.

  

**Options**

•  --batch_size: Number of SMILES per batch (max 100). Default is 100.

•  --max_retries: Maximum number of retries for failed batches. Default is 3.

•  --retry_delay: Delay between retries in seconds. Default is 10.

**Example**

```bash
	pyclassyfire data/unique_smiles.tsv results/ --batch_size 50 --max_retries 5 --retry_delay 15
```

## Handling Interruptions

If the classification process is disrupted or interrupted, you can easily resume by rerunning the same command with the original input file and output directory. The script will automatically identify any missing SMILES and continue processing from where it left off.

### Steps to Resume

1. **Rerun the Classification Command:**

    ```bash
    pyclassyfire classify data/unique_smiles.tsv results/ --batch_size 50 --max_retries 5 --retry_delay 15
    ```

2. **Automatic Detection:**

    The script will detect unprocessed SMILES and handle them accordingly, ensuring that previously completed batches are not reprocessed.

### Tips

- **Adjusting Batch Size:**

    If you encounter server overload or performance issues, consider reducing the `--batch_size` parameter to decrease the number of SMILES processed per batch. For example:

    ```bash
    pyclassyfire classify data/unique_smiles.tsv results/ --batch_size 25 --max_retries 5 --retry_delay 15
    ```

## Output Directory Structure

When you run the PyClassyFire classification script, the specified output directory will be organized as follows:

- **logs/**  
  This subdirectory contains log files that record the details of the classification process.


- **intermediate_results/**  
  Stores intermediate JSON files with results from each batch processed by the ClassyFire API.  Allows the script to resume from the last successful batch in case of interruptions.


- **output.json**  
  The final consolidated JSON file that merges all intermediate results. This file contains classification results for all processed SMILES strings, as returned by the ClassyFire API.


- **missing_smiles.json**  
  This file is generated to identify any problematic SMILES strings that were not successfully processed or matched to SMILES of ClassyFireAPI.  
  
- **Structure of `missing_smiles.json`:**
  
  ```json
  {
      "original_smiles_1": "canonical_smiles_1",
      "original_smiles_2": "canonical_smiles_2",
      ...
  }
  
For more detailed instructions and tutorials, refer to the [Notebooks](notebooks/) folder.

## Input File Format

The input file should be a TSV (Tab-Separated Values) file with a single column containing SMILES strings. Here’s an example of how the input file should look:

```tsv
SMILES
CCO
C1=CC=CC=C1
CC(C)C[C@@H](C(=O)O)NC(=O)N1CCC2=CC(=C(C=C2C1)OC)OC
...
```


•  **Header:** The first line should be a header, e.g., SMILES.

•  **SMILES Strings:** Each subsequent line contains a SMILES string representing a chemical compound.

**SMILES Canonicalization**

PyClassyFire canonicalizes SMILES strings using the RDKit library to ensure consistency. The canonization process utilizes the following RDKit function:
```python
Chem.MolToSmiles(self.mol, canonical=True)
```

**Output JSON Format**

The output JSON file contains the classification results as returned by the ClassyFire API. Each entry includes the original SMILES string and its corresponding classification information. SMILES in the output are formatted as per the ClassyFire API’s response structure.

Example of an output entry:
```json
[
    {
        "identifier": "Q12021439-1",
        "smiles": "COC1=CC=C(C=C1)C1=CN2C(=NC3=C2C(=O)N(C)C(=O)N3C)N1C1=C2OCCOC2=CC=C1",
        "inchikey": "InChIKey=NZQDSIFUFRNVJQ-UHFFFAOYSA-N",
        "kingdom": {
            "name": "Organic compounds",
            ...
     }
 ]
```

## Tutorials

  

For comprehensive tutorials and examples on using PyClassyFire, refer to the [Notebooks](notebooks/) folder. These Jupyter notebooks provide step-by-step guidance on setting up, processing data, and interpreting results.

## Acknowledgements

  

PyClassyFire was inspired by and builds upon the work of the following GitHub repositories:

•  [JamesJeffryes/pyclassyfire](https://github.com/JamesJeffryes/pyclassyfire.git)

•  [wykswr/classyfire_cli](https://github.com/wykswr/classyfire_cli.git)

  

We thank the authors for their valuable contributions and inspiration.

# To-Do List

- [ ] Task 1 (problematic smiles, unify canonization between api and ours, affect only small portion)
- [ ] Task 2 (verify installation)

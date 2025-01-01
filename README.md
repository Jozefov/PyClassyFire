# Classyfire CLI
This is a command line interface for the Classyfire API, which is based on the batch API calling. 
It allows you to classify chemical compounds at large scale.

## Installation
```bash
pip install classyfire-cli
```

## Usage
Use SMILES:
```bash
classyfire input.txt output.json
```

Use InChI:
```bash
classyfire --identifier InChI input.txt output.json
```

Export to TSV:
```bash
classyfire input.txt output.tsv
```
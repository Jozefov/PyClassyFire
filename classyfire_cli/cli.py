import click
import json
from .src.batch import Job, Scheduler
from .src.utils import MoleCule, chunk_tasks

@click.command()
@click.option('--identifier', type=click.Choice(['SMILES', 'InChI']), help='Type of identifier', default='SMILES')
@click.argument('input', type=click.File('r'))
@click.argument('output', type=click.File('w'))
def main(identifier, input, output):
    """
    Classify molecules using ClassyFire API
    """
    # get suffix for the input file
    input_suffix = input.name.split('.')[-1]
    if input_suffix == 'json':
        identifiers = json.load(input)
    elif input_suffix == 'txt':
        identifiers = input.readlines()
    else:
        raise ValueError('Input file must be a json or txt file')
    identifiers = [x.strip() for x in identifiers]
    if identifier == 'SMILES':
        smiles_list = [MoleCule.from_smiles(x).canonical_smiles for x in identifiers]
    else:
        smiles_list = [MoleCule.from_inchi(x).canonical_smiles for x in identifiers]
    if len(smiles_list) > 100:
        jobs = [Job(chunk) for chunk in chunk_tasks(smiles_list, 100)]
    else:
        jobs = [Job(smiles_list)]
    scheduler = Scheduler(jobs)
    click.echo('Running jobs...')
    scheduler.run()
    click.echo('Exporting results...')
    result = scheduler.export()
    smiles2identifier = dict(zip(smiles_list, identifiers))
    result = {smiles2identifier.get(k, k): v for k, v in result.items()}
    output_suffix = output.name.split('.')[-1]
    if output_suffix == 'json':
        json.dump(result, output)
    elif output_suffix == 'tsv':
        output.write('Identifier\tSuperclass\tClass\tSubclass\n')
        for k, v in result.items():
            output.write(f'{k}\t{v["superclass"]}\t{v["class"]}\t{v["subclass"]}\n')
    else:
        raise ValueError('Output file must be a json or tsv file')
    click.echo('Done')


if __name__ == '__main__':
    main()
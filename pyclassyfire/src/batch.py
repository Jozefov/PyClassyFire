from .api import get_results, structure_query
from .utils import take_class, MoleCule, load_existing_results, save_intermediate_results, chunk_tasks
import json
import time
import click
from tqdm import tqdm
from requests.exceptions import HTTPError, ConnectionError
from http.client import RemoteDisconnected
import logging
from datetime import datetime
import os

class Job:
    def __init__(self, smiles: list[str]):
        assert len(smiles) <= 100
        self.smiles = smiles
        self.query_id = None
        self.start_time = None

    def submit(self):
        # Use actual newline characters instead of literal backslashes
        self.query_id = structure_query("\n".join(self.smiles))
        self.start_time = time.time()

    @property
    def is_done(self) -> bool:
        result = json.loads(get_results(self.query_id))
        return result["classification_status"] == "Done"

    @property
    def is_stale(self) -> bool:
        return time.time() - self.start_time > 60 * 10  # 10 minutes

    def parse_results(self) -> list[dict]:
        result = json.loads(get_results(self.query_id))
        return result["entities"]


class Scheduler:
    """
    Only one job can be submitted at a time. When a job is running, other jobs are queued.
    """

    def __init__(self, jobs: list[Job]):
        self.jobs = jobs
        self.results = {}
        self.retry = 0

    def run(self):
        pbar = tqdm(total=len(self.jobs))
        while self.jobs:
            try:
                job = self.jobs[0]
                if job.query_id is None:
                    job.submit()
                    time.sleep(60)
                elif job.is_done:
                    time.sleep(10)
                    self.results[job.query_id] = job.parse_results()
                    self.jobs.pop(0)
                    pbar.update(1)
                    time.sleep(10)
                elif job.is_stale:
                    self.results[job.query_id] = []
                    self.jobs.pop(0)
                    pbar.update(1)
                else:
                    time.sleep(60)
            except HTTPError:
                if self.retry < 3:
                    self.retry += 1
                    time.sleep(300)
                else:
                    self.jobs.pop(0)
                    self.results[job.query_id] = []
                    pbar.update(1)
                    self.retry = 0

        pbar.close()

    def _aggregate(self) -> list[dict]:
        aggregation = []
        for _, results in self.results.items():
            aggregation.extend(results)
        return aggregation

    def export(self) -> dict[str, dict]:
        result = {}
        for record in self._aggregate():
            smiles = MoleCule.from_smiles(record['smiles']).canonical_smiles
            superclasses = take_class(record['superclass'])
            classes = take_class(record['class'])
            subclasses = take_class(record['subclass'])
            result[smiles] = {
                "superclass": superclasses,
                "class": classes,
                "subclass": subclasses
            }
        return result


def process_batches_with_saving_and_retry(
        smiles_list,
        batch_size=100,
        output_dir='../data/intermediate_results/',
        max_retries=3,
        retry_delay=10
):
    """
    Processes SMILES in batches with resumption and improved error handling.
    Saves each batch immediately after processing, including original SMILES for order-based matching.

    Parameters:
    - smiles_list (list): List of canonical SMILES strings to process.
    - batch_size (int): Number of SMILES per batch (max 100).
    - output_dir (str): Directory to save intermediate JSON files.
    - max_retries (int): Maximum number of retries for failed batches.
    - retry_delay (int): Delay between retries in seconds.

    Returns:
    - list: List of paths to saved intermediate JSON files.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Setup logging
    logs_dir = os.path.join(output_dir, 'logs')
    os.makedirs(logs_dir, exist_ok=True)

    # Generate unique log filename based on current time
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"log_{current_time}.log"
    log_filepath = os.path.join(logs_dir, log_filename)

    # Configure logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Prevent adding multiple handlers
    if not logger.handlers:
        fh = logging.FileHandler(log_filepath)
        fh.setLevel(logging.INFO)

        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)

        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)

        logger.addHandler(fh)
        logger.addHandler(ch)

    smiles_list = list(set(smiles_list))
    print(f'All SMILES: {len(smiles_list)}')
    logger.info(f'All SMILES: {len(smiles_list)}')

    # Load already processed results and the highest batch number
    already_processed, max_batch_num = load_existing_results(output_dir)
    print(f'Already processed SMILES: {len(already_processed)}')
    logger.info(f'Already processed SMILES: {len(already_processed)}')

    # Normalize SMILES returned from ClassyFire
    already_processed_tmp = []
    for smi in already_processed:
        canonical_smi = MoleCule.from_smiles(smi).canonical_smiles
        already_processed_tmp.append(canonical_smi)
    already_processed_set = set(already_processed_tmp)

    # Identify remaining SMILES to process
    remaining_smiles = set(smiles_list) - already_processed_set
    print(f'Remaining SMILES to process: {len(remaining_smiles)}')
    logger.info(f'Remaining SMILES to process: {len(remaining_smiles)}')

    # Remove duplicates from remaining_smiles
    remaining_smiles = list(set(remaining_smiles))
    print(f'Remaining unique SMILES to process after removing duplicates: {len(remaining_smiles)}')
    logger.info(f'Remaining unique SMILES to process after removing duplicates: {len(remaining_smiles)}')

    # Enforce batch size limit
    if batch_size > 100:
        print("Batch size cannot exceed 100. Setting batch_size to 100.")
        logger.warning("Batch size cannot exceed 100. Setting batch_size to 100.")
        batch_size = 100

    # Create batches
    batches = chunk_tasks(remaining_smiles, batch_size)
    total_batches = len(batches)

    print(f"Total remaining batches to process: {total_batches}")
    logger.info(f"Total remaining batches to process: {total_batches}")

    if total_batches == 0:
        print("All batches have already been processed.")
        return []

    saved_files = []

    pbar = tqdm(total=total_batches, desc="Processing Batches")

    for batch in batches:
        max_batch_num += 1
        batch_id = max_batch_num
        retries = 0
        while retries <= max_retries:
            try:
                # Initialize and submit the job
                job = Job(batch)
                job.submit()
                logging.info(f"Submitted Batch {batch_id} with Query ID {job.query_id}")
                print(f"Submitted Batch {batch_id} with Query ID {job.query_id}")

                # Poll for job status
                while True:
                    time.sleep(30)  # Wait before checking status
                    try:
                        result = get_results(job.query_id)
                        if not result:
                            raise ValueError("Empty response from API.")
                        result_json = json.loads(result)
                    except json.JSONDecodeError as e:
                        logging.error(f"JSON decode error for Batch {batch_id}: {e}")
                        raise e

                    classification_status = result_json.get("classification_status")
                    if classification_status == "Done":
                        molecules = result_json.get("entities", [])
                        expected_count = len(batch)
                        returned_count = len(molecules)

                        if returned_count != expected_count:
                            raise ValueError(f"Batch {batch_id}: Expected {expected_count} molecules, but received {returned_count}.")

                        if not molecules:
                            print(f"No results returned for Batch {batch_id}.")
                            logging.warning(f"No results returned for Batch {batch_id}.")
                        else:
                            print(f"Batch {batch_id} completed with {returned_count} molecules.")
                            logging.info(f"Batch {batch_id} completed with {returned_count} molecules.")

                        # Save intermediate results with original SMILES
                        save_intermediate_results(batch_id, molecules, batch, output_dir)
                        saved_files.append(os.path.join(output_dir, f'intermediate_{batch_id}.json'))
                        pbar.update(1)
                        break
                    elif classification_status == "Error":
                        error_message = result_json.get("error_message", "Unknown error.")
                        raise Exception(f"API Error for Batch {batch_id}: {error_message}")
                    else:
                        print(f"Batch {batch_id} status: {classification_status}. Waiting...")
                        logging.info(f"Batch {batch_id} status: {classification_status}. Waiting...")

                # If completed successfully, break out of the retry loop
                break

            except (click.exceptions.BadParameter, ConnectionError, Exception) as e:
                retries += 1
                if retries > max_retries:
                    print(f"Batch {batch_id}: Maximum retries reached. Skipping batch.")
                    logging.error(f"Batch {batch_id}: Maximum retries reached. Error: {e}")
                    # Save intermediate results with empty molecules and original SMILES
                    save_intermediate_results(batch_id, [], batch, output_dir)
                    saved_files.append(os.path.join(output_dir, f'intermediate_{batch_id}.json'))
                    pbar.update(1)
                    break
                else:
                    print(
                        f"Batch {batch_id}: Error encountered: {e}. Retrying ({retries}/{max_retries}) after {retry_delay} seconds...")
                    logging.error(
                        f"Batch {batch_id}: Error encountered: {e}. Retrying ({retries}/{max_retries}) after {retry_delay} seconds...")
                    time.sleep(retry_delay)

    pbar.close()
    logger.info("Processing completed.")
    print("Processing completed.")
    return saved_files
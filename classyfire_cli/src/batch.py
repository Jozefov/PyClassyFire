from .api import get_results, structure_query
from .utils import take_class, MoleCule
import json
import time
from tqdm import tqdm
from requests.exceptions import HTTPError


class Job:
    def __init__(self, smiles: list[str]):
        assert len(smiles) <= 1000
        self.smiles = smiles
        self.query_id = None
        self.start_time = None

    def submit(self):
        self.query_id = structure_query("\\n".join(self.smiles))
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

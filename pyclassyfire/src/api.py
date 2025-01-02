# api.py
import requests
import json

url = "http://classyfire.wishartlab.com"

def structure_query(compound, label='pyclassyfire'):
    """
    Submit a compound information to the ClassyFire service for evaluation
    and receive an ID which can be used to collect results.

    :param compound: The compound structures as line-delimited InChIKey or SMILES.
    :type compound: str
    :param label: A label for the query
    :type label: str
    :return: A query ID number
    :rtype: int
    """
    payload = {
        "label": label,
        "query_input": compound,
        "query_type": "STRUCTURE"
    }
    headers = {"Content-Type": "application/json"}
    r = requests.post(f"{url}/queries.json", data=json.dumps(payload), headers=headers)
    r.raise_for_status()
    return r.json()['id']

def get_results(query_id, return_format="json"):
    """
    Given a query_id, fetch the classification results.

    :param query_id: A numeric query ID returned at the time of query submission
    :type query_id: str
    :param return_format: Desired return format. Valid types are json, csv, or sdf
    :type return_format: str
    :return: Query information
    :rtype: str
    """
    r = requests.get(f"{url}/queries/{query_id}.{return_format}", headers={"Content-Type": f"application/{return_format}"})
    r.raise_for_status()
    return r.text
import tests.utilities as tu


"""
Normalize Variant Utility Tests
--------------------------------
"""


def test_normalize_variant_utility(client):
    url = tu.normalize_variant_query('variant=NM_021960.4%3Ac.740C>T')
    response = client.get(url)

    assert response.status_code == 200


def test_normalize_variant_utility_transcript_ref_disagree_position_variant(client):
    """
    This is a specific corner case where liftover fails if HGVS strict bounds checks are enforced
    Details here: https://github.com/biocommons/hgvs/issues/717
    """
    url = tu.normalize_variant_query('variant=NC_000001.10:g.145592073A>T')

    response = client.get(url)

    assert response.status_code == 200

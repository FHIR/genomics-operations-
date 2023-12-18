import tests.utilities as tu


"""
Normalize Variant Utility Tests
--------------------------------
"""


def test_normalize_variant_utility(client):
    url = tu.normalize_variant_query('variant=NM_021960.4%3Ac.740C>T')
    response = client.get(url)

    assert response.status_code == 200

import tests.utilities as tu


"""
Find Population Specific Variants Tests
---------------------------------------
"""


def test_population_specific_variants_variants_supplied(client):
    url = tu.find_population_specific_variants_query('genomicSourceClass=germline')
    response = client.get(url)

    assert response.status_code == 400


def test_population_specific_variants_1(client):
    url = tu.find_population_specific_variants_query('variants=NC_000011.10:g.8263343T>C&variants=NC_000006.11:26091305:T:C')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_SPECIFIC_VARIANTS_OUTPUT_DIR}1.json', response.json)


def test_population_specific_variants_2(client):
    # PDE4DIP gene is on minus strand in build 37, but is on positive strand in build 38.
    # The hgvs library returns "HGVSDataNotAvailableError: No alignments for NM_001002811.1 in GRCh38 using splign" in
    # this case.

    url = tu.find_population_specific_variants_query(
        'subject=HG00403&variants=NC_000001.10:144931726:G:A&variants=NC_000001.10:145532548:T:C&variants=NC_000001.10:145592072:A:T&includePatientList=true')
    response = client.get(url)

    assert response.status_code == 422


# def test_population_specific_variants_3(client):
#     url = tu.find_population_specific_variants_query('variants=NC_000001.10:144931726:G:A, NC_000001.10:145532548:T:C, NC_000001.10:145592072:A:T')
#     response = client.get(url)

#     tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_SPECIFIC_VARIANTS_OUTPUT_DIR}3.json', response.json)


"""
Find Population Structural Intersecting Variants Tests
------------------------------------------------------
"""


def test_population_structural_intersecting_variants_ranges_supplied(client):
    url = tu.find_population_structural_intersecting_variants_query('genomicSourceClass=germline')
    response = client.get(url)

    assert response.status_code == 400


def test_population_structural_intersecting_variants_1(client):
    url = tu.find_population_structural_intersecting_variants_query('ranges=NC_000022.10:42522500-42526812&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_STRUCTURAL_INTERSECTING_VARIANTS_OUTPUT_DIR}1.json', response.json)


def test_population_structural_intersecting_variants_2(client):
    url = tu.find_population_structural_intersecting_variants_query('ranges=NC_000017.11:43044294-43125364,NC_000013.11:32315507-32400268&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_STRUCTURAL_INTERSECTING_VARIANTS_OUTPUT_DIR}2.json', response.json)


def test_population_structural_intersecting_variants_3(client):
    url = tu.find_population_structural_intersecting_variants_query('ranges=NC_000022.10:42522500-42526812&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_STRUCTURAL_INTERSECTING_VARIANTS_OUTPUT_DIR}3.json', response.json)


"""
Find Population Structural Subsuming Variants Tests
---------------------------------------------------
"""


def test_population_structural_subsuming_variants_ranges_supplied(client):
    url = tu.find_population_structural_subsuming_variants_query('genomicSourceClass=germline')
    response = client.get(url)

    assert response.status_code == 400


def test_population_structural_subsuming_variants_1(client):
    url = tu.find_population_structural_subsuming_variants_query('ranges=NC_000022.10:42522500-42526812&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_STRUCTURAL_SUBSUMING_VARIANTS_OUTPUT_DIR}1.json', response.json)


def test_population_structural_subsuming_variants_2(client):
    url = tu.find_population_structural_subsuming_variants_query('ranges=NC_000017.10:37844346-37884911,NC_000007.13:116312249-116438431&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_STRUCTURAL_SUBSUMING_VARIANTS_OUTPUT_DIR}2.json', response.json)


"""
Find Population Specific Haplotypes Tests
-----------------------------------------
"""


def test_population_specific_haplotypes_haplotypes_supplied(client):
    url = tu.find_population_specific_haplotypes_query('genomicSourceClass=germline')
    response = client.get(url)

    assert response.status_code == 400


def test_population_specific_haplotypes_1(client):
    url = tu.find_population_specific_haplotypes_query('haplotypes=CYP2C9 *1/*2&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_SPECIFIC_HAPLOTYPES_OUTPUT_DIR}1.json', response.json)


def test_population_specific_haplotypes_2(client):
    url = tu.find_population_specific_haplotypes_query('haplotypes=CYP2C19&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_SPECIFIC_HAPLOTYPES_OUTPUT_DIR}2.json', response.json)


def test_population_specific_haplotypes_3(client):
    url = tu.find_population_specific_haplotypes_query('haplotypes=HLA-A*01&includePatientList=true&haplotypes=HLA-B*08&haplotypes=HLA-DRB1*03')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_SPECIFIC_HAPLOTYPES_OUTPUT_DIR}3.json', response.json)


def test_population_specific_haplotypes_4(client):
    url = tu.find_population_specific_haplotypes_query('haplotypes=HLA-A*01,HLA-B*08,HLA-DRB1*03&includePatientList=true')
    response = client.get(url)

    tu.compare_actual_and_expected_output(f'{tu.FIND_POPULATION_SPECIFIC_HAPLOTYPES_OUTPUT_DIR}4.json', response.json)

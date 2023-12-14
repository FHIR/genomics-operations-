# HL7 FHIR Genomics Operations - Reference Implementation
Source code for a [public reference implementation](https://fhir-gen-ops.herokuapp.com/) of [HL7 FHIR Genomics Operations](http://build.fhir.org/ig/HL7/genomics-reporting/operations.html).

Please refer to [project Wiki page](https://github.com/FHIR/genomics-operations/wiki) for details of this reference implementation, including how to replicate.

For additional information on the operations and the reference implementation, please see our [JAMIA manuscript](https://academic.oup.com/jamia/advance-article/doi/10.1093/jamia/ocac246/6957062).

Issues (bugs, enhancements, etc) can be entered [here](https://github.com/FHIR/genomics-operations/issues). (Legacy issues are [here](https://docs.google.com/spreadsheets/d/1xPRDB2lvMPTImPHLwUvSboILZLG6jH1LHVXoOfgak9U/edit#gid=0)). Contact info@elimu.io for other comments.

## Use Case

A common use case driving the operations is the notion of an application (e.g. a SMART-ON-FHIR clinical genomics App, a clinical decision support application, an EHR screen) needing specific genotype or phenotype information, for a patient or a population. Applications have diverse needs, such as matching a cancer patient to available clinical trials based on identified somatic variants; screening for actionable hereditary conditions; identifying a risk for adverse medication reactions based on pharmacogenomic variants; updating a patient's risk as knowledge of their variants evolves; and more. A goal for FHIR Genomics operations is to ultimately support any and all of these clinical scenarios.

## Scope

In scope are clinical genomics operations. In the future, operations supporting variant calling and annotation, and knowledge base lookups may be added. We further categorize clinical genomics operations along two orthogonal axes - subject vs. population, and genotype vs. phenotype. For example, the '**find-subject-variants**' operation retrieves genotype information for a single subject; whereas the '**find-population-tx-implications**' retrieves a count or list of patients having specific phenotypes (such as being intermediate metabolizers of clopidogrel).

## Response

All operations return a JSON output. However, if an invalid request is submitted, or some other error occurs, a JSON response is returned in the following format:

```javascript
{
  "type": string,
  "title": string,
  "detail": string,
  "status": int
}
```

## Status Codes

The operations return the following status codes:

| Status Code | Description |
| :--- | :--- |
| 200 | `Successfully executed request` |
| 400 | `ERROR: Invalid query parameters` |
| 404 | `ERROR: Patient not found` |
| 422 | `ERROR: Failed LiftOver` |
| 500 | `INTERNAL SERVER ERROR` |

## Testing

To run the [integration tests](https://github.com/FHIR/genomics-operations/tree/main/tests), you can use the VS Code Testing functionality which should discover them automatically. You can also
run `python3 -m pytest` from the terminal to execute them all.

Additionally, since the tests run against the Mongo DB database, if you need to update the test data in this repo, you
can run `OVERWRITE_TEST_EXPECTED_DATA=true python3 -m pytest` from the terminal and then create a pull request with the
changes.

## Heroku Deployment

Currently, there are two environments running in Heroku:
- Dev: <https://fhir-gen-ops-dev-ca42373833b6.herokuapp.com/>
- Prod: <https://fhir-gen-ops.herokuapp.com/>

Pull requests will trigger a deployment to the dev environment automatically after being merged.

The ["Manual Deployment"](https://github.com/FHIR/genomics-operations/actions/workflows/manual_deployment.yml) workflow
can be used to deploy code to either the `dev` or `prod` environments. To do so, please select "Run workflow", ignore
the "Use workflow from" dropdown which lists the branches in the current repo (I can't disable / remove it) and then
select the environment, the branch and the repository. By default, the `https://github.com/FHIR/genomics-operations`
repo is specified, but you can replace it with any any fork.

Deployments to the prod environment can only be triggered manually from the `main` branch of the repo using the Manual
Deployment.

### Heroku Stack

Make sure that the Python version under [`runtime.txt`](./runtime.txt) is
[supported](https://devcenter.heroku.com/articles/python-support#supported-runtimes) by the
[Heroku stack](https://devcenter.heroku.com/articles/stack) that is currently running in each environment.

### UTA Database

The Biocommons [hgvs](https://github.com/biocommons/hgvs) library which is used for variant parsing, validation and
normalisation requires access to a copy of the [UTA](https://github.com/biocommons/uta) Postgres database.

We have provisioned a Heroku Postgres instance in the Prod environment which contains the imported data from a database
dump as described [here](https://github.com/biocommons/uta#installing-from-database-dumps).

The connection string for this database can be found in Heroku under the `UTA_DATABASE_URL` environment variable.

Additionally, we define a `UTA_DATABASE_SCHEMA` environment variable in the [`.env`](.env) file which contains the name
of the currently imported database schema.

Database import procedure (it will take about 10 minutes):

```shell
> UTA_SCHEMA="uta_20210129b" # Specify the UTA schema you wish to use
> PGPASSWORD="${POSTGRES_PASSWORD}"
> gzip -cdq ${UTA_SCHEMA}.pgd.gz | grep -v anonymous | psql -U ${POSTGRES_USER} -1 -v ON_ERROR_STOP=1 -d ${POSTGRES_DATABASE} -h ${POSTGRES_HOST} -Eae
```

Note: `grep -v anonymous` is required because it's not possible to create an `anonymous` role in Heroku Postgres.

Once the process finishes, if you are using the Heroku Postgres Basic plan on the
[Essential Tier](https://devcenter.heroku.com/articles/heroku-postgres-plans#essential-tier), you'll bump into the 10
million rows / database limit. However, it's safe to ignore the warnings about this limit, since Heroku will simply
revoke INSERT privileges from the database and the hgvs library only needs read-only access to this database.

### RefSeq access

The RefSeq metadata from the UTA database needs to be in sync with the RefSeq data which is available for the *Seqfetcher
Utility* endpoint. See the [`fetchRefseq.py` utility docs](./utilities/fetchRefseq.py) for details.

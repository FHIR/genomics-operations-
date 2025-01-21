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

For local development, you will have to create a `secrets.env` file in the root of the repo and add in it the MongoDB
password and, optionally, the UTA Postgres database connection string (see the UTA section below for details):

```
MONGODB_READONLY_PASSWORD=...
UTA_DATABASE_URL=...
```

To run the [integration tests](https://github.com/FHIR/genomics-operations/tree/main/tests), you can use the VS Code
Testing functionality which should discover them automatically. You can also run `python3 -m pytest` from the terminal
to execute them all.

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

We define a `UTA_DATABASE_SCHEMA` environment variable in the [`.env`](.env) file which contains the name of the
currently imported database schema.

#### Database import procedure (it will take about 30 minutes):

- Go to the UTA dump download site (http://dl.biocommons.org/uta/) and get the latest `<UTA_SCHEMA>.pgd.gz` file.
- Go to https://dashboard.heroku.com/apps/fhir-gen-ops/resources and click on the "Heroku Postgres" instance (it will
open a new window)
- Go to the Settings tab
- Click "View Credentials"
- Use the fields from this window to fill in the variables below

```shell
$ POSTGRES_HOST="<Heroku Postgres Host>"
$ POSTGRES_DATABASE="<Heroku Postgres Database>"
$ POSTGRES_USER="<Heroku Postgres User>"
$ PGPASSWORD="<Heroku Postgres Password>"
$ UTA_SCHEMA="<UTA Schema>" # Specify the UTA schema of the UTA dump you downloaded (example: uta_20240523b)
$ gzip -cdq ${UTA_SCHEMA}.pgd.gz | grep -v '^GRANT USAGE ON SCHEMA .* TO anonymous;$' | grep -v '^ALTER .* OWNER TO uta_admin;$' | psql -U ${POSTGRES_USER} -1 -v ON_ERROR_STOP=1 -d ${POSTGRES_DATABASE} -h ${POSTGRES_HOST} -Eae
```

Note: The `grep -v` commands are required because the Heroku Postgres instance doesn't allow us to create a new role.

Once complete, make sure you update the `UTA_DATABASE_SCHEMA` environment variable in the [`.env`](.env) file and commit
it.

#### Connection string

The connection string for this database can be found in the same Heroku Postgres Settings tab under "View Credentials".
It is pre-populated in the Heroku runtime under the `UTA_DATABASE_URL` environment variable. Additionally, we set the
same `UTA_DATABASE_URL` environment variable in GitHub so the CI can can use this database when running the tests.

For local development, if you'd like to use this Postgres instance instead of the HGVS public one
(`postgresql://anonymous:anonymous@uta.biocommons.org/uta`), please add `UTA_DATABASE_URL` with the Heroku Postgres
connection string in the `secrets.env` file.

#### Testing the database

```shell
$ pgcli "${UTA_DATABASE_URL}"
> set schema '<UTA Schema>'; # Specify the UTA schema of the UTA dump you downloaded (example: uta_20240523b)
> select count(*) from alembic_version
    union select count(*) from associated_accessions
    union select count(*) from exon
    union select count(*) from exon_aln
    union select count(*) from exon_set
    union select count(*) from gene
    union select count(*) from meta
    union select count(*) from origin
    union select count(*) from seq
    union select count(*) from seq_anno
    union select count(*) from transcript
    union select count(*) from translation_exception;
```

### RefSeq data

The RefSeq metadata from the UTA database needs to be in sync with the RefSeq data which is available for the *Seqfetcher
Utility* endpoint. Currently, this is stored in GitHub as release artifacts.

To update the RefSeq data, you will have to install `seqrepo` locally and run `./utilities/pack_seqrepo_data.py`. Here
is a step-by-step guide on how to do this:

```shell
$ mkdir seqrepo
$ cd seqrepo
$ python3 -m venv .venv
$ . .venv/bin/activate
$ pip install setuptools==75.7.0
$ pip install biocommons.seqrepo==0.6.9
$ # See https://github.com/biocommons/biocommons.seqrepo/issues/171 for a bug that's causing issues with the builtin
$ # rsync on OSX.
$ brew install rsync # OSX-specific. Guess the standard package managers have it available on Linux
$ seqrepo --rsync-exe /opt/homebrew/bin/rsync -r . pull --update-latest
$ # If you'll get a "Permission denied" error, then you can run the following command (using the temp directory which got created):
$ # > chmod +w 2024-02-20.r4521u5y && mv 2024-02-20.r4521u5y 2024-02-20 && ln -s 2024-02-20 latest
$
$ # cd to genomics-operations repo
$ python ./utilities/pack_seqrepo_data.py --seqrepo_dir /path/to/seqrepo/dir/latest
$ # Upload tar archives from ./tmp/ to a new GitHub release and then update `UTILITIES_DATA_VERSION` in the `.env` file
$ # such that it contains the short SHA of the new release which contains the updated data.
```

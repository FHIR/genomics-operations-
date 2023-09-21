from pathlib import Path
from urllib.request import urlretrieve
import tarfile
import sys

SEQ_FILES_RELEASE_TAG = '113c119'


def download_refseq_files():
    try:
        # Make sure the parent folder exists
        Path('refseq').mkdir(exist_ok=True)

        for build in ['GRCh37', 'GRCh38']:
            filename = build + 'seq.tar.gz'
            filepath = f'refseq/{filename}'

            if not Path(filepath).is_file():
                # Download seq files
                urlretrieve('https://github.com/FHIR/genomics-operations/releases/download/' + SEQ_FILES_RELEASE_TAG + '/' +
                            filename, filepath)

            # Unarchive seq files each time on startup, just in case something went wrong the previous time
            with tarfile.open(filepath, "r|gz") as f:
                f.extractall(path='refseq')

    except Exception as err:
        sys.exit(f"Failed to download refseq files: {err=}, {type(err)=}")

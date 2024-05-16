import connexion
import flask
from flask_cors import CORS
import os

import hgvs
# Disable the hgvs LRU cache to avoid blowing up memory
# TODO: Revisit this, since this caching might not use a ton of memory.
hgvs.global_config.lru_cache.maxsize = 0
# Disable HGVS strict bounds checks as a workaround for liftover failures: https://github.com/biocommons/hgvs/issues/717
hgvs.global_config.mapping.strict_bounds = False


def create_app():
    # App and API
    options = {
        'swagger_url': '/',
        'openapi_spec_path': '/openapi.json'
    }

    connex_app = connexion.App("__name__", specification_dir='app/', options=options)

    api = connex_app.add_api('api_spec.yml')

    bp = flask.Blueprint('redoc', __name__, url_prefix=api.base_path,
                         template_folder=os.path.dirname(__file__))

    specPath = api.base_path + api.options.openapi_spec_path
    def serveRedoc(): return flask.render_template('redoc.j2', openapi_spec_url=specPath)

    bp.add_url_rule('/docs/', __name__, serveRedoc)

    app = connex_app.app
    app.register_blueprint(bp)

    # Allow CORS for all domains on all routes
    CORS(app, resources={
        r"/subject-operations/.*": {"origins": "*"},
        r"/population-operations/.*": {"origins": "*"},
        r"/utilities/.*": {"origins": "*"}
    })

    # Avoid sorting the output JSON
    app.config['JSON_SORT_KEYS'] = False

    return app


if __name__ == '__main__':
    create_app().run()

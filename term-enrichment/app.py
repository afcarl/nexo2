from term_enrichment import TermEnrichment

from flask_restful import Api, Resource
from flask import Flask


class TopResource(Resource):

    def get(self):
        message = {
            'message': 'NeXO Term Enrichment server.  Use any of the following endpoints.',
            'endpoints': [
                '/enrich'
            ]
        }
        return message


if __name__ == '__main__':

    app = Flask(__name__)
    api = Api(app)

    # Endpoints
    api.add_resource(TopResource, '/')
    api.add_resource(TermEnrichment, '/enrich')

    app.run(host='0.0.0.0', debug=True)

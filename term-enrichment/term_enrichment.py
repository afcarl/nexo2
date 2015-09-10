from flask_restful import Resource, reqparse
from hypergeometric import HypergeometricTest

GENES = 'genes'
ALPHA = 'alpha'
ONTOLOGY_TYPE = 'type'

# Valid namespace
VALID_ONTOLOGY = ['NEXO', 'MF', 'CC', 'BP']


class TermEnrichment(Resource):

    def __init__(self):
        self.__parser = reqparse.RequestParser()
        self.__parser.add_argument(GENES, type=str)
        self.__parser.add_argument(ALPHA, type=float)
        self.__parser.add_argument(ONTOLOGY_TYPE, type=str)

        self.__tester = HypergeometricTest()

    def get(self):
        return {'results': 'use POST to get results.'}

    def post(self):
        results = []
        args = self.__parser.parse_args()

        # Gene names are case INSENSITIVE!
        genes = args[GENES]
        cutoff = args[ALPHA]
        ontology_type = args[ONTOLOGY_TYPE]

        if ontology_type is None:
            ontology_type = 'NEXO'
        elif ontology_type not in VALID_ONTOLOGY:
            return {'results': 'Ontology type not supported: ' + ontology_type}

        gene_list = genes.split()
        genes_upper = []
        for gene in gene_list:
            genes_upper.append(gene.upper())

        # print(genes_upper)
        return self.__tester.perform_test(ontology_type, genes_upper, cutoff)

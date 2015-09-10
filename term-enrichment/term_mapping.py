# root terms
ROOT_NEXO = "joining_root"
ROOT_BP = "GO:0008150"
ROOT_MF = "GO:0003674"
ROOT_CC = "GO:0005575"


class TermMapping:

    _gene2terms = {}
    _term2genes = {}
    _num_genes = {}

    def __init__(self):
        print "Loading resource files..."

        # NeXO resources
        self._gene2terms['NEXO'] = self.__create_map("data/gene2terms.txt", ",")
        self._term2genes['NEXO'] = self.__create_map("data/term2genes.txt", ",")

        # GO resources
        self._gene2terms['BP'] = self.__create_map("data/biological_process.info_gain.gene_term.genes", "|")
        self._term2genes['BP'] = self.__create_map("data/biological_process.info_gain.gene_term.terms", "|")

        self._gene2terms['CC'] = self.__create_map("data/cellular_component.info_gain.gene_term.genes", "|")
        self._term2genes['CC'] = self.__create_map("data/cellular_component.info_gain.gene_term.terms", "|")

        self._gene2terms['MF'] = self.__create_map("data/molecular_function.info_gain.gene_term.genes", "|")
        self._term2genes['MF'] = self.__create_map("data/molecular_function.info_gain.gene_term.terms", "|")

        # Count genes.  Root term has all genes.
        all_genes = self._term2genes['NEXO'][ROOT_NEXO]
        self._num_genes['NEXO'] = len(all_genes)

        all_genes = self._term2genes['MF'][ROOT_MF]
        self._num_genes['MF'] = len(all_genes)
        all_genes = self._term2genes['CC'][ROOT_CC]
        self._num_genes['CC'] = len(all_genes)
        all_genes = self._term2genes['BP'][ROOT_BP]
        self._num_genes['BP'] = len(all_genes)

        print("Ready.  NeXO Genes = " + str(self._num_genes['NEXO']))
        print("MF Genes = " + str(self._num_genes['MF']))
        print("BP Genes = " + str(self._num_genes['BP']))
        print("CC Genes = " + str(self._num_genes['CC']))

    def __create_map(self, file_name, delimiter):
        mapping = {}
        with open(file_name, "r") as f:
            for line in f:
                line = line.rstrip()
                parts = line.split("\t")
                gene_id = parts[0]
                terms = parts[1].split(delimiter)
                mapping[gene_id] = terms
        f.close()
        return mapping

    def get_gene_mapping(self, ontology_type):
        return self._gene2terms[ontology_type]

    def get_term_mapping(self, ontology_type):
        return self._term2genes[ontology_type]

    def get_gene_count(self, ontology_type):
        return self._num_genes[ontology_type]

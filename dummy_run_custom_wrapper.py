import confinspect
from confinspect import NetworkInferenceWrapper, TestRunner, Selectors, InputHandler
import os

class CustomWrapper(NetworkInferenceWrapper.NetworkInferenceWrapper):

    def _infer_network(self, expression_data, rank):
        print('infer network from a single block!')
        pass

    def _get_top_k_edges(self, i, k):
        print('Get top k edges, with each edge being a tuple of nodes, and don\'t forget to order the nodes alphabetically, if the \
            edges are undirected!')
        pass

if __name__ == "__main__":
    InputHandler.setup_directories()
    wrap = CustomWrapper(os.getcwd())
    data = {'HNSC': {'tcga': True, 'ged': f'TCGA-HNSC.htseq_fpkm.tsv', 'pt': f'TCGA-HNSC.GDC_phenotype.tsv', 'sep': ',', 'tissue_type_field': None, 'tissue_type': None}}
    fields = {'gender.demographic': {'role': Selectors.Role.CONFOUNDER, 'type': Selectors.BlockType.CATEGORY}}
    params = {'N_from': 0, 'N_to': 1, 'M_from': 0, 'M_to': 1, 'k_max': 100, 'save_networks': False, 'combine': False, 'g_all': False, 'logfile': 'log.txt'}
    InputHandler.verify_input(data, params, fields)
    tr = TestRunner.TestRunner(data, fields, params)
    tr.add_custom_algorithm(wrap, 'name')
    tr.induce_partitions()
    print(tr.rnd_partitions['HNSC']['gender.demographic'][0]['male'])
    print(len(tr.rnd_partitions['HNSC']['gender.demographic'][0]['male']))

    print(set(tr.rnd_partitions['HNSC']['gender.demographic'][0]['male']).intersection(set(tr.rnd_partitions['HNSC']['gender.demographic'][0]['female'])))

    print('This will throw an error unless the CustomWrapper gets implemented properly.')
    tr.run_all()


class BioThingsHelper:
    def get_gene(gene_id):
        gene=gene_client.getgene(gene_id, fields='symbol,name')
        return gene
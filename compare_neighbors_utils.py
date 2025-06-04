def compare_gene_neighbors(up_a, down_a, up_b, down_b):
    genes_a = [g["gene"] for g in up_a + down_a]
    genes_b = [g["gene"] for g in up_b + down_b]

    set_a = set(genes_a)
    set_b = set(genes_b)

    shared = set_a.intersection(set_b)
    unique_a = set_a - set_b
    unique_b = set_b - set_a

    return sorted(shared), sorted(unique_a), sorted(unique_b)
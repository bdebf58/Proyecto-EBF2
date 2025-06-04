import gzip

def extract_neighbors(gff_file, gene_id, flank_genes=3, is_gz=False):
    opener = gzip.open if is_gz else lambda f, mode: f
    features = []

    with opener(gff_file, 'rt') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attributes = fields
            if feature == "gene":
                gene_name = ""
                for attr in attributes.split(";"):
                    if attr.startswith("ID="):
                        gene_name = attr.replace("ID=", "").split(":")[-1]
                        break
                features.append({
                    "gene": gene_name,
                    "seqid": seqid,
                    "start": int(start),
                    "end": int(end),
                    "attributes": attributes
                })

    # Ordenar por coordenada
    features.sort(key=lambda x: x["start"])

    # Buscar el índice del ortólogo
    index = next((i for i, feat in enumerate(features) if gene_id in feat["gene"]), None)
    if index is None:
        return None, None, None

    neighbors_upstream = features[max(0, index - flank_genes): index]
    neighbors_downstream = features[index + 1: index + 1 + flank_genes]

    return features[index], neighbors_upstream, neighbors_downstream
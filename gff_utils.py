import gzip

def extract_exons_utr3(file_obj, gene_id, is_gz=False):
    exons = []
    utr3s = []
    opener = gzip.open if is_gz else lambda f, mode: f

    with opener(file_obj, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            if gene_id in attributes:
                if feature_type == "exon":
                    exons.append((int(start), int(end)))
                elif feature_type == "three_prime_UTR":
                    utr3s.append((int(start), int(end)))

    num_exons = len(exons)
    utr3_length = sum(end - start + 1 for start, end in utr3s)
    return num_exons, utr3_length
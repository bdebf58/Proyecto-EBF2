import gzip
from io import BytesIO, TextIOWrapper

def is_gzipped_by_content(file_bytes):
    return file_bytes[:2] == b'\x1f\x8b'

def find_gene_id_from_protein(uploaded_file, protein_id):
    if uploaded_file is None:
        return None

    file_bytes = uploaded_file.read()
    buffer = BytesIO(file_bytes)

    if is_gzipped_by_content(file_bytes):
        f = gzip.open(buffer, 'rt')
    else:
        f = TextIOWrapper(buffer, encoding="utf-8")

    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        attributes = fields[8]
        if protein_id in attributes:
            for part in attributes.split(";"):
                if part.startswith("Parent="):
                    return part.replace("Parent=", "").split(":")[-1]
    return None
import streamlit as st
import sqlite3
import pandas as pd
import plotly.express as px
import os
from BCBio import GFF
from Bio import SeqIO
import gzip
import io
import re

DB_PATH = "species_records.db"

# Conexión y helpers
def get_connection():
    return sqlite3.connect(DB_PATH)

def get_all_records():
    conn = get_connection()
    df = pd.read_sql_query("SELECT * FROM Records", conn)
    conn.close()
    return df

def update_record(record_id, field, value):
    conn = get_connection()
    c = conn.cursor()
    c.execute(f"UPDATE Records SET {field} = ? WHERE id = ?", (value, record_id))
    conn.commit()
    conn.close()

def export_csv(df):
    return df.to_csv(index=False).encode('utf-8')

def listar_archivos_desde_data(extensiones):
    carpeta = "data"
    if not os.path.exists(carpeta):
        return []
    return [f for f in os.listdir(carpeta) if any(f.endswith(ext) for ext in extensiones)]

# Helper para cargar archivo (de streamlit u open local)
def cargar_archivo_con_opcion(label, extensiones, key_prefix):
    fuente = st.radio(f"Selecciona origen del archivo para {label}:", ["Subir desde PC", "Elegir desde carpeta 'data'"], key=key_prefix+"_fuente")

    archivo = None
    archivo_nombre = None

    if fuente == "Subir desde PC":
        archivo = st.file_uploader(f"Sube archivo {label} aquí", type=extensiones, key=key_prefix+"_upload")
        if archivo is not None:
            archivo_nombre = archivo.name
    else:
        archivos = listar_archivos_desde_data(extensiones)
        archivo_sel = st.selectbox(f"Selecciona archivo {label} desde carpeta 'data'", archivos, key=key_prefix+"_select")
        if archivo_sel:
            ruta = os.path.join("data", archivo_sel)
            archivo = open(ruta, "rb")
            archivo_nombre = archivo_sel

    return archivo, archivo_nombre, fuente


# Interfaz Streamlit
st.set_page_config(layout="wide")
st.title("🧬Base de Datos EBF2")

tab1, tab2, tab3 = st.tabs(["📋 Ver tabla", "✏️ Editar registros", "📊 Visualización & Exportación"])

# TABLA COMPLETA
with tab1:
    st.subheader("📋 Todos los registros")
    df = get_all_records()
    st.dataframe(df, use_container_width=True)

# EDICIÓN DE REGISTROS
with tab2:
    st.subheader("✏️ Editar un registro existente")
    df = get_all_records()
    record_ids = df['id'].tolist()
    selected_id = st.selectbox("Selecciona ID del registro a editar", record_ids)
    selected_row = df[df['id'] == selected_id].iloc[0]

    field = st.selectbox("Campo a modificar", df.columns.drop('id'))
    new_value = st.text_input("Nuevo valor", str(selected_row[field]))
    if st.button("Actualizar"):
        update_record(selected_id, field, new_value)
        st.success(f"Registro {selected_id} actualizado correctamente.")

# VISUALIZACIONES Y EXPORTACIÓN
with tab3:
    st.subheader("📊 Visualización de datos")
    df = get_all_records()

    col1, col2 = st.columns(2)
    with col1:
        fig1 = px.histogram(df, x="clade", title="Distribución por Clado")
        st.plotly_chart(fig1, use_container_width=True)

    with col2:
        if "percent_identity" in df.columns:
            try:
                df['percent_identity'] = pd.to_numeric(df['percent_identity'], errors='coerce')
                fig2 = px.box(df.dropna(subset=["percent_identity"]), x="clade", y="percent_identity", title="% Identidad por Clado")
                st.plotly_chart(fig2, use_container_width=True)
            except Exception as e:
                st.warning("No se pudo generar gráfico de identidad.")

    st.subheader("📥 Exportar registros")
    st.download_button("📤 Descargar CSV", export_csv(df), file_name="EBF2_records_export.csv", mime="text/csv")



from gff_utils import extract_exons_utr3

with st.expander("📥 Analizar exones y 3'UTR desde archivo GFF3 (.gff3/.gz)"):
    gff_file, gff_file_name, fuente = cargar_archivo_con_opcion("GFF3 (.gff3/.gz)", [ "gff3", "gz"], "gff_exons")
    gene_id_input = st.text_input("ID del gen ortólogo (ej: AT5G25350)", key="gene_input")

    if gff_file and gene_id_input:
        is_gz = gff_file_name.endswith(".gz") if gff_file_name else False

        num_exons, utr3_len = extract_exons_utr3(gff_file, gene_id_input, is_gz=is_gz)
        st.success(f"Número de exones: {num_exons} | Longitud del 3'UTR: {utr3_len} nt")

        if st.button("📌 Guardar en base de datos"):
            conn = get_connection()
            c = conn.cursor()
            c.execute("UPDATE Records SET exons = ?, utr3 = ? WHERE ortholog_id = ?", (num_exons, utr3_len, gene_id_input))
            conn.commit()
            conn.close()
            st.success("Base de datos actualizada con éxito.")

            if fuente == "Elegir desde carpeta 'data'":
                gff_file.close()  # cerrar archivo abierto localmente



from synteny_utils import extract_neighbors

with st.expander("🧬 Analizar genes vecinos (sintenia)"):
    gff_file_synt, gff_file_synt_name, fuente_synt = cargar_archivo_con_opcion("GFF3 (.gff3/.gz) para sintenia", ["gff3", "gz"], "gff_synt")
    gene_id_synt = st.text_input("ID del gen ortólogo a comparar (ej: AT5G25350)", key="gene_synt")

    if gff_file_synt and gene_id_synt:
        is_gz = gff_file_synt_name.endswith(".gz") if gff_file_synt_name else False
        center, upstream, downstream = extract_neighbors(gff_file_synt, gene_id_synt, is_gz=is_gz)

        if center:
            st.info(f"🧬 Gen central: {center['gene']} ({center['start']} - {center['end']})")
            st.markdown("**⬆️ Genes upstream:**")
            for u in upstream:
                st.write(f"↖️ {u['gene']} ({u['start']} - {u['end']})")
            st.markdown("**⬇️ Genes downstream:**")
            for d in downstream:
                st.write(f"↘️ {d['gene']} ({d['start']} - {d['end']})")
        else:
            st.error("No se encontró el gen en el archivo GFF.")

        if fuente_synt == "Elegir desde carpeta 'data'":
            gff_file_synt.close()


from protein_to_gene_utils import find_gene_id_from_protein

with st.expander("🔍 Buscar ID de gen desde ID de proteína (archivo GFF)"):
    gff_protein_file, gff_protein_file_name, fuente_prot = cargar_archivo_con_opcion("GFF3 (.gff3/.gz) para búsqueda", ["gff3", "gz"], "gff_protein")
    protein_id_input = st.text_input("ID de proteína ortóloga (ej: XP_006842065.1)", key="prot_input")

    if gff_protein_file and protein_id_input:
        gene_id_result = find_gene_id_from_protein(gff_protein_file, protein_id_input)
        if gene_id_result:
            st.success(f"Gene ID correspondiente: {gene_id_result}")
        else:
            st.error("No se encontró un gene ID correspondiente en el archivo.")

        if fuente_prot == "Elegir desde carpeta 'data'":
            gff_protein_file.close()

with st.expander("🔄 Buscar gene ID y analizar genes vecinos automáticamente"):
    gff_combo_file, gff_combo_file_name, fuente_combo = cargar_archivo_con_opcion("GFF3 (.gff3/.gz)", ["gff3", "gz"], "gff_combo")
    protein_input = st.text_input("ID de proteína ortóloga (ej: XP_006842065.1)", key="prot_combo_input")

    if gff_combo_file and protein_input:
        found_gene = find_gene_id_from_protein(gff_combo_file, protein_input)

        if found_gene:
            st.success(f"Gene ID correspondiente: {found_gene}")
            is_gz = gff_combo_file_name.endswith(".gz") if gff_combo_file_name else False
            center, upstream, downstream = extract_neighbors(gff_combo_file, found_gene, is_gz=is_gz)

            if center:
                st.info(f"🧬 Gen central: {center['gene']} ({center['start']} - {center['end']})")
                st.markdown("**⬆️ Genes upstream:**")
                for u in upstream:
                    st.write(f"↖️ {u['gene']} ({u['start']} - {u['end']})")
                st.markdown("**⬇️ Genes downstream:**")
                for d in downstream:
                    st.write(f"↘️ {d['gene']} ({d['start']} - {d['end']})")
            else:
                st.warning("Gene ID encontrado, pero no se pudo analizar vecinos.")
        else:
            st.error("No se encontró un gene ID para ese ID de proteína.")

        if fuente_combo == "Elegir desde carpeta 'data'":
            gff_combo_file.close()


with st.expander("🧬 Comparar genes vecinos entre dos especies (sintenia comparada)"):
    col1, col2 = st.columns(2)
    with col1:
        gff_species_a, gff_species_a_name, fuente_a = cargar_archivo_con_opcion("GFF3 especie A (referencia)", ["gff3", "gz"], "gff_species_a")
        gene_a = st.text_input("🧬 ID de gen en especie A (ej: AT5G25350)", key="gene_a")
    with col2:
        gff_species_b, gff_species_b_name, fuente_b = cargar_archivo_con_opcion("GFF3 especie B (ortólogo)", ["gff3", "gz"], "gff_species_b")
        gene_b = st.text_input("🧬 ID de gen ortólogo en especie B", key="gene_b")

    if gff_species_a and gff_species_b and gene_a and gene_b:
        is_gz_a = gff_species_a_name.endswith(".gz") if gff_species_a_name else False
        is_gz_b = gff_species_b_name.endswith(".gz") if gff_species_b_name else False

        center_a, upstream_a, downstream_a = extract_neighbors(gff_species_a, gene_a, is_gz=is_gz_a)
        center_b, upstream_b, downstream_b = extract_neighbors(gff_species_b, gene_b, is_gz=is_gz_b)

        if center_a and center_b:
            genes_a = [g["gene"] for g in upstream_a + downstream_a]
            genes_b = [g["gene"] for g in upstream_b + downstream_b]
            set_a = set(genes_a)
            set_b = set(genes_b)
            shared = set_a.intersection(set_b)
            unique_a = set_a - set_b
            unique_b = set_b - set_a

            st.success(f"🔗 Genes vecinos compartidos: {len(shared)}")
            st.markdown("**🟢 Comunes:**")
            for g in sorted(shared):
                st.write(f"✔️ {g}")

            st.markdown("**🔴 Únicos en especie A:**")
            for g in sorted(unique_a):
                st.write(f"❌ {g}")

            st.markdown("**🔵 Únicos en especie B:**")
            for g in sorted(unique_b):
                st.write(f"❌ {g}")
        else:
            st.error("No se pudieron obtener los genes vecinos de una de las especies.")

        if fuente_a == "Elegir desde carpeta 'data'":
            gff_species_a.close()
        if fuente_b == "Elegir desde carpeta 'data'":
            gff_species_b.close()


from compare_neighbors_utils import compare_gene_neighbors
from visual_synteny import plot_synteny_tracks

with st.expander("🧬 Comparar genes vecinos + visualización"):
    col1, col2 = st.columns(2)
    with col1:
        gff_species_a, gff_species_a_name, fuente_a = cargar_archivo_con_opcion("GFF3 especie A (referencia)", ["gff3", "gz"], "gff_va")
        gene_a = st.text_input("ID gen especie A (ej: AT5G25350)", key="gene_viz_a")
    with col2:
        gff_species_b, gff_species_b_name, fuente_b = cargar_archivo_con_opcion("GFF3 especie B (ortólogo)", ["gff3", "gz"], "gff_viz_b")
        gene_b = st.text_input("ID gen especie B (ej: Eucgr.Hxxxx)", key="gene_viz_b")

    if gff_species_a and gff_species_b and gene_a and gene_b:
        is_gz_a = gff_species_a_name.endswith(".gz") if gff_species_a_name else False
        is_gz_b = gff_species_b_name.endswith(".gz") if gff_species_b_name else False

        dfa, dfb, df_union = compare_gene_neighbors(gff_species_a, gff_species_b, gene_a, gene_b, is_gz_a, is_gz_b)
        fig = plot_synteny_tracks(dfa, dfb, df_union)
        st.plotly_chart(fig, use_container_width=True)

        if fuente_a == "Elegir desde carpeta 'data'":
            gff_species_a.close()
        if fuente_b == "Elegir desde carpeta 'data'":
            gff_species_b.close()


def cargar_fasta(file, is_gz=True):
    if is_gz:
        return SeqIO.to_dict(SeqIO.parse(gzip.open(file, "rt"), "fasta"))
    else:
        return SeqIO.to_dict(SeqIO.parse(file, "fasta"))

def parse_gff3(gff_file, gene_id, is_gz=True):
    exons = []
    chrom, gene_start, gene_end, strand = None, None, None, "+"
    found_gene = False

    open_func = gzip.open if is_gz else open
    with open_func(gff_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue
            seqid, source, feature_type, start, end, score, strand_, phase, attributes = parts

            attr_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in attributes.split(";") if "=" in kv}
            if not found_gene:
                if (feature_type == "gene" and ("ID" in attr_dict and gene_id in attr_dict["ID"])):
                    chrom = seqid
                    gene_start = int(start)
                    gene_end = int(end)
                    strand = strand_
                    found_gene = True
            if found_gene and feature_type == "exon":
                parent = attr_dict.get("Parent", "")
                if gene_id in parent:
                    exons.append((int(start), int(end)))

    return chrom, gene_start, gene_end, strand, sorted(exons)

def get_sequence(chrom, start, end, strand, fasta_dict):
    seq = fasta_dict.get(chrom, None)
    if not seq:
        return None
    region = seq.seq[start - 1:end]
    return str(region.reverse_complement()) if strand == "-" else str(region)

def get_protein_sequence(gene_id, protein_dict):
    for record in protein_dict.values():
        if gene_id in record.description or gene_id in record.id:
            return str(record.seq)
    return None

def extraer_secuencia(gff_file, genome_fasta, protein_fasta, gene_id, tipo, is_gz=True):
    fasta_dict = cargar_fasta(genome_fasta, is_gz)
    protein_dict = cargar_fasta(protein_fasta, is_gz)

    chrom, start, end, strand, exons = parse_gff3(gff_file, gene_id, is_gz)
    if chrom is None:
        return None

    if tipo == "DNA":
        return get_sequence(chrom, start, end, strand, fasta_dict)
    elif tipo == "mRNA":
        mseq = ""
        for ex_start, ex_end in exons:
            mseq += get_sequence(chrom, ex_start, ex_end, strand, fasta_dict)
        return mseq
    elif tipo == "Proteína":
        return get_protein_sequence(gene_id, protein_dict)

    return None


with st.expander("🧬 Extraer secuencias del gen"):
    gff_file = st.file_uploader("📄 Archivo GFF3", type=["gff3", "gz"])
    genome_fasta = st.file_uploader("🧬 FASTA de genoma (DNA)", type=["fa", "fasta", "gz"], key="genome")
    protein_fasta = st.file_uploader("🧪 FASTA de proteínas", type=["fa", "fasta", "gz"], key="protein")
    gene_id = st.text_input("🔍 ID del gen (ej: AT5G25350)")
    tipo = st.selectbox("🧬 Tipo de secuencia a extraer", ["DNA", "mRNA", "Proteína"])

    if gff_file and genome_fasta and protein_fasta and gene_id:
        seq = extraer_secuencia(gff_file, genome_fasta, protein_fasta, gene_id, tipo, is_gz=True)
        if seq:
            st.code(seq, language="text")
        else:
            st.error("No se pudo encontrar o extraer la secuencia.")

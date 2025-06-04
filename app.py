import streamlit as st
import sqlite3
import pandas as pd
import plotly.express as px

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

# Interfaz Streamlit
st.set_page_config(layout="wide")
st.title("🧬 Proyecto Rafael - Base de Datos EBF2")

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
    gff_file = st.file_uploader("Sube archivo .gff3 o .gff3.gz", type=["gff3", "gz"], key="gff_exons")
    gene_id_input = st.text_input("ID del gen ortólogo (ej: AT5G25350)", key="gene_input")

    if gff_file and gene_id_input:
        is_gz = gff_file.name.endswith(".gz")
        num_exons, utr3_len = extract_exons_utr3(gff_file, gene_id_input, is_gz=is_gz)
        st.success(f"Número de exones: {num_exons} | Longitud del 3'UTR: {utr3_len} nt")

        if st.button("📌 Guardar en base de datos"):
            conn = get_connection()
            c = conn.cursor()
            c.execute("UPDATE Records SET exons = ?, utr3 = ? WHERE ortholog_id = ?", (num_exons, utr3_len, gene_id_input))
            conn.commit()
            conn.close()
            st.success("Base de datos actualizada con éxito.")


from synteny_utils import extract_neighbors

with st.expander("🧬 Analizar genes vecinos (sintenia)"):

    gff_file_synt = st.file_uploader("Sube archivo .gff3 o .gff3.gz para sintenia", type=["gff3", "gz"], key="gff_synt")
    gene_id_synt = st.text_input("ID del gen ortólogo a comparar (ej: AT5G25350)", key="gene_synt")

    if gff_file_synt and gene_id_synt:
        is_gz = gff_file_synt.name.endswith(".gz")
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


from protein_to_gene_utils import find_gene_id_from_protein

with st.expander("🔍 Buscar ID de gen desde ID de proteína (archivo GFF)"):
    gff_protein_file = st.file_uploader("Sube archivo .gff3 o .gz para búsqueda", type=["gff3", "gz"], key="gff_protein")
    protein_id_input = st.text_input("ID de proteína ortóloga (ej: XP_006842065.1)", key="prot_input")

    if gff_protein_file and protein_id_input:
        is_gz = gff_protein_file.name.endswith(".gz")
        if gff_protein_file is not None and protein_id_input:
            gene_id_result = find_gene_id_from_protein(gff_protein_file, protein_id_input)
            if gene_id_result:
                st.success(f"Gene ID correspondiente: {gene_id_result}")
            else:
                st.error("No se encontró un gene ID correspondiente en el archivo.")
            st.success(f"Gene ID correspondiente: {gene_id_result}")
        else:
            st.error("No se encontró un gene ID correspondiente en el archivo.")


with st.expander("🔄 Buscar gene ID y analizar genes vecinos automáticamente"):
    gff_combo_file = st.file_uploader("Sube archivo .gff3 o .gz", type=["gff3", "gz"], key="gff_combo")
    protein_input = st.text_input("ID de proteína ortóloga (ej: XP_006842065.1)", key="prot_combo_input")

if gff_combo_file and protein_input:
    found_gene = find_gene_id_from_protein(gff_combo_file, protein_input)

    if found_gene:
        st.success(f"Gene ID correspondiente: {found_gene}")
        is_gz = gff_combo_file.name.endswith(".gz")
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



with st.expander("🧬 Comparar genes vecinos entre dos especies (sintenia comparada)"):
    col1, col2 = st.columns(2)
    with col1:
        gff_species_a = st.file_uploader("📂 GFF3 especie A (referencia, ej. Arabidopsis)", type=["gff3", "gz"], key="gff_species_a")
        gene_a = st.text_input("🧬 ID de gen en especie A (ej: AT5G25350)", key="gene_a")
    with col2:
        gff_species_b = st.file_uploader("📂 GFF3 especie B (ortólogo)", type=["gff3", "gz"], key="gff_species_b")
        gene_b = st.text_input("🧬 ID de gen ortólogo en especie B", key="gene_b")

    if gff_species_a and gff_species_b and gene_a and gene_b:
        is_gz_a = gff_species_a.name.endswith(".gz")
        is_gz_b = gff_species_b.name.endswith(".gz")

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

from compare_neighbors_utils import compare_gene_neighbors


from visual_synteny import plot_synteny_tracks

with st.expander("🧬 Comparar genes vecinos + visualización"):
    col1, col2 = st.columns(2)
    with col1:
        gff_species_a = st.file_uploader("📂 GFF3 especie A (referencia)", type=["gff3", "gz"], key="gff_va")
        gene_a = st.text_input("ID gen especie A", key="gva")
    with col2:
        gff_species_b = st.file_uploader("📂 GFF3 especie B (ortóloga)", type=["gff3", "gz"], key="gff_vb")
        gene_b = st.text_input("ID gen especie B", key="gvb")

    if gff_species_a and gff_species_b and gene_a and gene_b:
        isgz_a = gff_species_a.name.endswith(".gz")
        isgz_b = gff_species_b.name.endswith(".gz")

        ca, upa, downa = extract_neighbors(gff_species_a, gene_a, is_gz=isgz_a)
        cb, upb, downb = extract_neighbors(gff_species_b, gene_b, is_gz=isgz_b)

        if ca and cb:
            shared, _, _ = compare_gene_neighbors(upa, downa, upb, downb)
            st.plotly_chart(plot_synteny_tracks(ca, upa, downa, cb, upb, downb, shared), use_container_width=True)
        else:
            st.warning("No se pudieron obtener vecinos de ambos genes.")

with st.expander("📂 Cargar archivos FASTA / GFF para analizar 3'UTR y exones"):
    uploaded_fasta = st.file_uploader("Sube archivo .fasta", type=["fasta", "fa"])
    uploaded_gff = st.file_uploader("Sube archivo .gff o .gff3", type=["gff", "gff3"])

    if uploaded_fasta and uploaded_gff:
        from io import StringIO
        from Bio import SeqIO

        gff_lines = uploaded_gff.read().decode("utf-8").splitlines()
        fasta_sequences = SeqIO.to_dict(SeqIO.parse(uploaded_fasta, "fasta"))

        utr_lengths = {}
        exon_counts = {}

        for line in gff_lines:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attributes = fields
            if feature == "three_prime_UTR":
                gene_id = attributes.split(";")[0].split("=")[-1]
                utr_len = int(end) - int(start) + 1
                utr_lengths[gene_id] = utr_lengths.get(gene_id, 0) + utr_len
            elif feature == "exon":
                gene_id = attributes.split(";")[0].split("=")[-1]
                exon_counts[gene_id] = exon_counts.get(gene_id, 0) + 1

        st.success(f"UTRs detectados: {len(utr_lengths)}, Exones detectados: {len(exon_counts)}")

        with st.expander("🔍 Mostrar resumen de genes analizados"):
            for gene_id in sorted(set(utr_lengths.keys()) | set(exon_counts.keys())):
                st.write(f"🧬 {gene_id} - 3'UTR: {utr_lengths.get(gene_id, 0)} nt, Exones: {exon_counts.get(gene_id, 0)}")

        if st.checkbox("Actualizar registros existentes con los datos cargados (si coinciden por ID)"):
            conn = get_connection()
            c = conn.cursor()
            for gene_id in utr_lengths:
                c.execute("UPDATE Records SET utr3 = ? WHERE ortholog_id = ?", (utr_lengths[gene_id], gene_id))
            for gene_id in exon_counts:
                c.execute("UPDATE Records SET exons = ? WHERE ortholog_id = ?", (exon_counts[gene_id], gene_id))
            conn.commit()
            conn.close()
            st.success("Registros actualizados con éxito desde los archivos proporcionados.")
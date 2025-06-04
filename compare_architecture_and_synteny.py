import streamlit as st
import pandas as pd
import plotly.express as px
from gff_utils import extract_exons_utr3
from synteny_utils import extract_neighbors
from compare_neighbors_utils import compare_gene_neighbors
from visual_synteny import plot_synteny_tracks
from app import get_all_records

# Punto 2: Comparación de arquitectura génica vs Arabidopsis
with st.expander("🧬 Comparar arquitectura génica con Arabidopsis"):
    df = get_all_records()
    ref_gene = df[df["ortholog_id"] == "AT5G25350"]
    if not ref_gene.empty:
        ref_exons = ref_gene.iloc[0]["exons"]
        ref_utr3 = ref_gene.iloc[0]["utr3"]

        df["exons_diff"] = df["exons"] - ref_exons
        df["utr3_diff"] = pd.to_numeric(df["utr3"], errors="coerce") - ref_utr3

        st.write("🧬 Diferencia en número de exones y longitud del 3'UTR respecto a Arabidopsis (EBF2)")
        st.dataframe(df[["species", "ortholog_id", "exons", "exons_diff", "utr3", "utr3_diff"]])

        fig = px.scatter(df, x="exons_diff", y="utr3_diff", hover_name="species", title="Comparación estructural vs Arabidopsis")
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.warning("No se encontró el gen AT5G25350 como referencia.")


# Punto 4: Comparación de sintenia con Arabidopsis
with st.expander("🧬 Ver si se conserva la sintenia con Arabidopsis"):
    arabidopsis_gff = st.file_uploader("📂 GFF de Arabidopsis", type=["gff3", "gz"], key="ref_synt")
    ortho_gff = st.file_uploader("📂 GFF de especie ortóloga", type=["gff3", "gz"], key="ortho_synt")
    ortho_gene_id = st.text_input("ID del gen ortólogo a comparar", key="ortho_gene_id")

    if arabidopsis_gff and ortho_gff and ortho_gene_id:
        ca, upa, downa = extract_neighbors(arabidopsis_gff, "AT5G25350", is_gz=arabidopsis_gff.name.endswith(".gz"))
        cb, upb, downb = extract_neighbors(ortho_gff, ortho_gene_id, is_gz=ortho_gff.name.endswith(".gz"))

        if ca and cb:
            shared, _, _ = compare_gene_neighbors(upa, downa, upb, downb)
            st.success(f"Genes vecinos compartidos: {len(shared)}")
            st.plotly_chart(plot_synteny_tracks(ca, upa, downa, cb, upb, downb, shared))
        else:
            st.error("No se pudo encontrar el gen o sus vecinos en alguno de los archivos.")
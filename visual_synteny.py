import plotly.graph_objects as go

def plot_synteny_tracks(center_a, up_a, down_a, center_b, up_b, down_b, shared_genes):
    def build_track(center, upstream, downstream, y, label_prefix):
        genes = upstream[::-1] + [center] + downstream
        x = 0
        positions = []
        for gene in genes:
            positions.append((gene["gene"], x))
            x += (gene["end"] - gene["start"]) // 1000 + 10
        return [
            go.Scatter(
                x=[p[1] for p in positions],
                y=[y]*len(positions),
                mode="markers+text",
                marker=dict(size=12, color=["green" if p[0] in shared_genes else "gray" for p in positions]),
                text=[p[0] for p in positions],
                textposition="bottom center",
                name=label_prefix
            )
        ]

    fig = go.Figure()
    fig.add_traces(build_track(center_a, up_a, down_a, 1, "Especie A"))
    fig.add_traces(build_track(center_b, up_b, down_b, 0, "Especie B"))

    fig.update_layout(
        title="🧬 Comparación visual de genes vecinos (sintenia)",
        yaxis=dict(
            tickvals=[1, 0],
            ticktext=["Especie A", "Especie B"]
        ),
        height=400,
        showlegend=False
    )
    return fig
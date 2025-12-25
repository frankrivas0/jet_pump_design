import plotly.graph_objects as go


def plot_data(results_table, ipr_table):
    fig = go.Figure()
    # PLOT OUTFLOW CURVE
    fig.add_trace(go.Scatter(
        y=results_table['Pressure@pwf'],
        x=results_table['Production bbl/d'],
        mode='lines+markers',
        name='Outflow Curve',
        line=dict(color='blue')
    ))
    # PLOT INFLOW CURVE
    fig.add_trace(go.Scatter(
        y=ipr_table['Inflow Bottomhole Pressure (psi)'],
        x=ipr_table['Inflow Rate (bbl/day)'],
        mode='lines+markers',
        name='Inflow Curve',
        line=dict(color='green')
    ))
    # PLOT CAVITATION CURVE
    fig.add_trace(go.Scatter(
        x=results_table['Cavitation flow bbl/d'],
        y=results_table['Pressure@pwf'],
        mode='lines+markers',
        fill='tozeroy',
        fillcolor='rgba(255,0,0,0.2)',
        name='Cavitation Curve',
        line=dict(color='red', dash='dash')
    ))

    fig.update_layout(
        title="Outflow and Inflow Curves",
        yaxis_title="Bottomhole Pressure (psi)",
        xaxis_title="Rate (bbl/day)",
        legend_title="Curves",
        template="plotly_white"
    )
    fig.update_yaxes(rangemode="nonnegative")

    return fig
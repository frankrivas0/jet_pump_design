import streamlit as st
from utils.generate_IPR_curve import ipr_curve
from utils.generate_pump_curve import pump_curve
from utils.plot import plot_data
from utils.find_operating_point import find_operating_point



st.set_page_config(page_title="Results", layout="wide")
if 'range_data_entered' not in st.session_state:
    st.session_state['range_data_entered'] = False

if 'well_data_uploaded' not in st.session_state:
    st.session_state['well_data_uploaded'] = False

if 'uploaded_file_name' in st.session_state:
    st.sidebar.text(f"Uploaded file: {st.session_state['uploaded_file_name']}")

# FORM TO INPUT BOTTOMHOLE PRESSURE RANGE FOR SIMULATION
with st.form("Bottomhole Pressure range"):
    st.header("Bottomhole Pressure Range for Simulation")
    pwfmin = st.number_input("Minimum Bottomhole Pressure (psi):", format="%.2f")
    pwfmax = st.number_input("Maximum Bottomhole Pressure (psi):", format="%.2f")
    range_data_entered = st.form_submit_button("Set Range")

# CHECK THAT THE RANGE IS VALID
if range_data_entered:
    if pwfmin <= 0 or pwfmax <= 0 or pwfmin >= pwfmax:
        st.error("Please enter a valid pressure range where min < max and both are greater than 0.")
        st.session_state['range_data_entered'] = False
        st.session_state['range_data'] = None
    else:
        st.success("Pressure range set!")
        st.session_state['range_data'] = {
            "pwfmin": pwfmin,
            "pwfmax": pwfmax
        }
        st.session_state['range_data_entered'] = True

simulate = st.button("Run Simulation")

if simulate:
    if st.session_state["well_data_uploaded"] and st.session_state["range_data_entered"]:
        outflow_results = pump_curve(st.session_state["well_data"], st.session_state["range_data"])
        inflow_results = ipr_curve(st.session_state["well_data"])
        st.session_state["outflow_results"] = outflow_results
        st.session_state["inflow_results"] = inflow_results
        fig = plot_data(outflow_results, inflow_results)
        st.session_state["plot_figure"] = fig
    else:
        st.error("No data available, incomplete or incorrect input.")

if "plot_figure" in st.session_state:
    st.header("Simulation Results")
    st.plotly_chart(st.session_state["plot_figure"], use_container_width=True)
    st.divider()
    with st.form("Pump Operating Point"):
        st.text("Based on the plot, enter a range in the x-axis where the outflow and inflow curves intersect:")
        x_low = st.number_input("Enter lower value:", format="%.2f")
        x_hi = st.number_input("Enter upper value:", format="%.2f")
        find_op = st.form_submit_button("Find Operating Point")
    if find_op:
        if x_low <= 0 or x_hi <= 0 or x_low >= x_hi:
            st.error("Please enter a valid range where min < max and both are greater than 0.")
        else:
            st.session_state['operating_point_dict'] = find_operating_point(st.session_state["outflow_results"], 
                                                                            st.session_state["inflow_results"], 
                                                                            x_low, 
                                                                            x_hi)
    if 'operating_point_dict' in st.session_state:
        st.subheader("Pump Operating Point Results")
        op_dict = st.session_state['operating_point_dict']
        st.write(f"**Rate (bbl/day):** {op_dict['Operating Rate (bbl/day)']:.2f}")
        st.write(f"**Bottomhole Pressure (psi):** {op_dict['Operating Bottomhole Pressure (psi)']:.2f}")
        st.write(f"**Efficiency (%):** {op_dict['Efficiency (%)']:.2f}")
        st.write(f"**Horsepower (HP):** {op_dict['Horsepower (HP)']:.2f}")

import streamlit as st

st.set_page_config(page_title="Table Results", layout="wide")

if 'uploaded_file_name' in st.session_state:
    st.sidebar.text(f"Uploaded file: {st.session_state['uploaded_file_name']}")
    
if 'outflow_results' in st.session_state and 'inflow_results' in st.session_state:
    st.header("Outflow Data Table")
    st.dataframe(st.session_state["outflow_results"], use_container_width=True)
    st.divider()
    st.header("Inflow Data Table")
    st.dataframe(st.session_state["inflow_results"], use_container_width=True)
else:
    st.warning("No simulation results available. Please run the simulation on the '📊 Plot Results' page.")
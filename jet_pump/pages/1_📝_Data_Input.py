import streamlit as st
from utils.schema import json_template, WellData
import json
from pydantic import ValidationError


st.set_page_config(page_title="Data Input", layout="wide")

if 'well_data_uploaded' not in st.session_state:
    st.session_state['well_data_uploaded'] = False

# PROVIDE A JSON TEMPLATE FOR THE USER TO FILL IN
st.sidebar.text("You can download a JSON template file to fill in your data or do it manually:")
st.sidebar.download_button(
    label="📥 Download JSON Data Template",
    data=json_template,
    file_name="data_template.json",
    mime="application/json"
)

# SECTION TO UPLOAD THE JSON FILE
st.sidebar.text("Please upload your filled JSON data file:")
uploaded_file = st.sidebar.file_uploader("Upload JSON Data File", type=["json"])

# CHECK IF A FILE HAS BEEN UPLOADED
if uploaded_file is not None:
    try:
        # Load JSON
        st.session_state['uploaded_file_name'] = uploaded_file.name
        data = json.load(uploaded_file)

        # Validate against schema
        well_data = WellData(**data)

        st.success("✅ JSON file is valid. You can proceed to the results page.")
        st.session_state['well_data'] = well_data.model_dump()
        st.session_state['well_data_uploaded'] = True

    except json.JSONDecodeError:
        st.session_state['well_data'] = None
        st.session_state['well_data_uploaded'] = False
        st.error("❌ Invalid JSON format")

    except ValidationError as e:
        st.session_state['well_data'] = None
        st.session_state['well_data_uploaded'] = False
        st.error("❌ JSON validation failed")
        st.json(e.errors())

if 'uploaded_file_name' in st.session_state:
    st.sidebar.text(f"Uploaded file: {st.session_state['uploaded_file_name']}")

# MANUAL DATA ENTRY SECTION
manual_entries = {}

with st.form("well_data_form"):
    col1, col2 = st.columns(2)
    with col1:
        st.header("Well Data")
        manual_entries["well_name"] = st.text_input("Well name:")
        manual_entries["perforations_depth_ft"] = st.number_input("Perforations depth (ft):", format="%.2f")
        manual_entries["pump_depth_ft"] = st.number_input("Pump depth (ft):", format="%.2f")
        manual_entries["tubingID_in"] = st.number_input("Tubing ID (in):", format="%.2f")
        manual_entries["tubingOD_in"] = st.number_input("Tubing OD (in):", format="%.2f")
        manual_entries["casingID_in"] = st.number_input("Casing ID (in):", format="%.2f")
        manual_entries["injection_pressure_psi"] = st.number_input("Injection pressure (psi):", format="%.2f")
        manual_entries["wellhead_pressure_psi"] = st.number_input("Wellhead pressure (psi):", format="%.2f")
        manual_entries["wellhead_temperature_degF"] = st.number_input("Wellhead temperature (°F):", format="%.2f")
        st.divider()

        st.header("Fluid Properties")
        manual_entries["api"] = st.number_input("Oil API gravity:", format="%.2f")
        manual_entries["bubble_point_pressure_psi"] = st.number_input("Oil Bubble point pressure (psi):", format="%.2f")
        manual_entries["gas_gravity"] = st.number_input("Gas specific gravity:", format="%.2f")
        manual_entries["water_gravity"] = st.number_input("Water specific gravity:", format="%.2f")
        manual_entries["injection_viscosity_cP"] = st.number_input("Injection fluid viscosity (cP):", format="%.2f")
    with col2:

        st.header("Production Data")
        manual_entries["production_rate_STBD"] = st.number_input("Production flow rate (STB/D):", format="%.2f")
        manual_entries["bsw"] = st.number_input("BSW (%):", format="%.2f")/100
        manual_entries["GOR"] = st.number_input("Gas-Oil Ratio (SCF/STB):", format="%.2f")
        st.divider()

        st.header("Reservoir Data")
        manual_entries["reservoir_pressure_psi"] = st.number_input("Reservoir pressure (psi):", format="%.2f")
        manual_entries["bottomhole_temperature_degF"] = st.number_input("Bottomhole temperature (°F):", format="%.2f")
        st.divider()

        st.header("Test Data")
        manual_entries["q_test"] = st.number_input("Test flow rate (STB/D):", format="%.2f")
        manual_entries["pwf_test"] = st.number_input("Test bottomhole pressure (psi):", format="%.2f")
        st.divider()

        st.header("Jet Pump Properties")
        manual_entries["pump_model"] = st.text_input("Pump name:")
        manual_entries["nozzle_area_in2"] = st.number_input("Nozzle area (in²):", format="%.5f")
        manual_entries["throat_area_in2"] = st.number_input("Throat area (in²):", format="%.5f")

    data = st.form_submit_button("Submit")

# SAVE ALL THE FIELD IMPUTS INTO A JSON OBJECT TO VALIDATE

if data:
    # CHECK THAT ALL FIELDS ARE FILLED AND EXPECTED TYPE
    try: 
        well_data = WellData(**manual_entries)
        st.success("✅ All fields are valid")
        st.session_state['well_data'] = well_data.model_dump()
        st.session_state['well_data_uploaded'] = True
    except ValidationError as e:
        st.session_state['well_data'] = None
        st.session_state['well_data_uploaded'] = False
        st.error("❌ Manual data validation failed")
        st.json(e.errors())

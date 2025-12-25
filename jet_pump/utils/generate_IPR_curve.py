from utils.jet import ipr
from pandas import DataFrame
from numpy import linspace, array

def ipr_curve(well_data, n_points=50):
    q_test = well_data["q_test"]
    pwf_test = well_data["pwf_test"]
    p_res = well_data["reservoir_pressure_psi"]
    pb = well_data["bubble_point_pressure_psi"]
    api = well_data["api"]
    h_perforations = well_data["perforations_depth_ft"]
    total_depth = well_data["pump_depth_ft"]

    x = ipr(q_test, pwf_test, p_res, pb)
    pwf = linspace(0, p_res, n_points)
    q_inflow = array([x.voguel(pwfi) for pwfi in pwf])
    pwf_corrected = pwf - 0.433*(141.5/(api+131.5))*(h_perforations - total_depth)
    inflow_curve = DataFrame({
        "Inflow Bottomhole Pressure (psi)": pwf_corrected,
        "Inflow Rate (bbl/day)": q_inflow
    })
    inflow_curve = inflow_curve[inflow_curve["Inflow Bottomhole Pressure (psi)"] >= 0]
    return inflow_curve
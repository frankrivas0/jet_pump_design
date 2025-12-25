from utils.jet import jet_pump
from pandas import DataFrame
from numpy import linspace, array

def pump_curve(well_data, range_data, n_points=50):
    p_inj = well_data["injection_pressure_psi"]
    p_wellhead = well_data["wellhead_pressure_psi"]
    pb = well_data["bubble_point_pressure_psi"]
    t_bottom = well_data["bottomhole_temperature_degF"]
    t_wellhead = well_data["wellhead_temperature_degF"]
    total_depth = well_data["pump_depth_ft"]
    dcsg = well_data["casingID_in"]
    dtbgID = well_data["tubingID_in"]
    dtbgOD = well_data["tubingOD_in"]
    api = well_data["api"]
    GOR = well_data["GOR"]
    bsw = well_data["bsw"]
    yg = well_data["gas_gravity"]
    yw = well_data["water_gravity"]
    mu_inj = well_data["injection_viscosity_cP"]
    q_test = well_data["q_test"]
    pwf_test = well_data["pwf_test"]
    aj = well_data["nozzle_area_in2"]
    at = well_data["throat_area_in2"]
    p_res = well_data["reservoir_pressure_psi"]
    pwfmin = range_data["pwfmin"]
    pwfmax = range_data["pwfmax"]
    
    pwf = linspace(pwfmin, pwfmax, n_points)

    q_prod = q_test  # Initial guess for production rate
    values = []
    for pwfi in pwf:
        result = jet_pump(pwfi, q_prod, aj, at, p_inj, p_wellhead, pb, t_bottom, 
                          t_wellhead, total_depth, dcsg, dtbgID, dtbgOD, api, GOR, bsw, yg, yw, mu_inj)
        values.append(result)

    results_table = DataFrame(values, columns = ['Production bbl/d', 'Pressure@pwf', 'Injection bbl/d', 
                                                 'Pressure@Discharge', 'Pressure@Nozzle', 'HP', 'Eff %', 
                                                 'Cavitation area sq-in', 'Cavitation flow bbl/d'])
    
    return results_table
        



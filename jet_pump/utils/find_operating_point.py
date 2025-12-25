from scipy.interpolate import interp1d

def find_operating_point(outflow_results, inflow_results, x_low, x_hi):
    # Extract relevant data
    q_inflow = inflow_results['Inflow Rate (bbl/day)']
    pwf_inflow = inflow_results['Inflow Bottomhole Pressure (psi)']
    q_outflow = outflow_results['Production bbl/d']
    pwf_outflow = outflow_results['Pressure@pwf']

    # Create interpolation functions
    f1 = interp1d(q_inflow, pwf_inflow, kind='cubic')
    f2 = interp1d(q_outflow, pwf_outflow, kind='linear')
    from scipy.optimize import brentq

    def difference(x):
        return f1(x) - f2(x)

    # Choose a range where you expect the intersection
    q_intersect = brentq(difference, x_low, x_hi)
    pwf_intersect = f1(q_intersect)

    # Find the efficiency at the operating point
    efficiency_interp = interp1d(outflow_results['Production bbl/d'], outflow_results['Eff %'], kind='linear')
    efficiency_at_op = efficiency_interp(q_intersect)

    # Find the horsepower at the operating point
    horsepower_interp = interp1d(outflow_results['Production bbl/d'], outflow_results['HP'], kind='linear')
    horsepower_at_op = horsepower_interp(q_intersect)
    return {
        "Operating Rate (bbl/day)": q_intersect,
        "Operating Bottomhole Pressure (psi)": pwf_intersect,
        "Efficiency (%)": efficiency_at_op,
        "Horsepower (HP)": horsepower_at_op
    }


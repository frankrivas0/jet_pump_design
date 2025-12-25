from pydantic import BaseModel, field_validator, Field

json_template = """{
    "well_name": "Test Well",
    "perforations_depth_ft": 9980,
    "pump_depth_ft": 9630,
    "tubingID_in": 2.992,
    "tubingOD_in": 3.5,
    "casingID_in": 6.276,
    "injection_pressure_psi": 3550,
    "wellhead_pressure_psi": 88,
    "wellhead_temperature_degF": 90,
    "bottomhole_temperature_degF": 205.0,
    "api": 29.1,
    "bubble_point_pressure_psi": 725,
    "gas_gravity": 1.4124,
    "water_gravity": 1.03,
    "injection_viscosity_cP": 2.5,
    "production_rate_STBD": 80,
    "bsw": 0.4,
    "GOR": 945,
    "reservoir_pressure_psi": 1151,
    "q_test": 134,
    "pwf_test": 530,
    "pump_model": "Sample Pump",
    "nozzle_area_in2": 0.0122,
    "throat_area_in2": 0.0311
}"""

class WellData(BaseModel):
    well_name: str = Field(..., description="Name of the well", examples=["Well A"])
    perforations_depth_ft: float = Field(..., gt=0, description="Depth of perforations in feet", examples=[5000.0])
    pump_depth_ft: float = Field(..., gt=0, description="Depth of the pump in feet", examples=[8000.0])
    tubingID_in: float = Field(..., gt=0, description="Internal diameter of tubing in inches", examples=[2.875])
    tubingOD_in: float = Field(..., gt=0, description="External diameter of tubing in inches", examples=[3.5])
    casingID_in: float = Field(..., gt=0, description="Internal diameter of casing in inches", examples=[5.5])
    injection_pressure_psi: float = Field(..., gt=0, description="Injection pressure in psi", examples=[2500.0])
    wellhead_pressure_psi: float = Field(..., gt=0, description="Wellhead pressure in psi", examples=[90.0])
    wellhead_temperature_degF: float = Field(..., description="Wellhead temperature in degrees Fahrenheit", examples=[150.0])
    bottomhole_temperature_degF: float = Field(..., description="Bottomhole temperature in degrees Fahrenheit", examples=[200.0])
    api: float = Field(..., gt=0, description="API gravity of the fluid", examples=[35.0])
    bubble_point_pressure_psi: float = Field(..., gt=0, description="Bubble point pressure in psi", examples=[2500.0])
    gas_gravity: float = Field(..., gt=0, description="Gas gravity", examples=[0.65])
    water_gravity: float = Field(..., gt=0, description="Water gravity", examples=[1.0])
    injection_viscosity_cP: float = Field(..., gt=0, description="Injection fluid viscosity in centipoise", examples=[1.0])
    production_rate_STBD: float = Field(..., gt=0, description="Production rate in STB/D", examples=[2000.0])
    bsw: float = Field(..., ge=0, le=1, description="Water cut", examples=[0.2])
    GOR: float = Field(..., ge=0, description="Gas-oil ratio", examples=[800.0])
    reservoir_pressure_psi: float = Field(..., gt=0, description="Reservoir pressure in psi", examples=[3000.0])
    q_test: float = Field(..., gt=0, description="Test production rate", examples=[1500.0])
    pwf_test: float = Field(..., gt=0, description="Test wellbore flowing pressure", examples=[1800.0])
    pump_model: str = Field(..., description="Model of the pump", examples=["Sample Pump"])
    nozzle_area_in2: float = Field(..., gt=0, description="Nozzle area in square inches", examples=[0.05])
    throat_area_in2: float = Field(..., gt=0, description="Throat area in square inches", examples=[0.1])

    class Config:
        extra = "forbid"

    @field_validator("well_name", "pump_model")
    def non_empty_string(cls, v):
        if not v or not v.strip():
            raise ValueError("Well name cannot be empty")
        return v
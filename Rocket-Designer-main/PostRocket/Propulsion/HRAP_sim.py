from hrap_python.api.core import sim_loop
from hrap_python.api.models import *
from hrap_python.api.sim import shift_OF
from hrap_python.util.units import *
import matplotlib.pyplot as plt


def TNAM(grain_diameter, port_diameter , grain_length, oxidizer_volume):

    material = MaterialData('hrap_python/propellants/Paraffin.npz')
    prop_config = PropellantConfig(ID=LengthValue(port_diameter, LengthUnit.METERS),
                                OD=LengthValue(grain_diameter, LengthUnit.METERS),
                                length=LengthValue(grain_length, LengthUnit.METERS),
                                cstar_eff=0.8,
                                material=material,
                                regression_model=shift_OF)

    hardware_config = HardwareConfig(tank_volume=VolumeValue(oxidizer_volume, VolumeUnit.CU_METERS),
                                    chamber_volume=None,
                                    injector_Cd=0.36,
                                    injector_D=LengthValue(0.375, LengthUnit.INCHES),
                                    injector_N=1,
                                    vent_state=VentConfig.NONE,
                                    vent_Cd=0.6,
                                    vent_D=LengthValue(0.028, LengthUnit.INCHES))

    nozzle_config = NozzleConfig(Cd=0.95,
                                throat=LengthValue(0.03, LengthUnit.METERS),
                                exp_ratio=1.88,
                                efficiency=0.95)

    sim_config = SimulationConfig(timestep=0.001,
                                max_run_time=150.0,
                                max_burn_time=15.0,
                                ambient_pressure=PressureValue(1.0, PressureUnit.ATM),
                                hardware=hardware_config,
                                prop=prop_config,
                                nozzle=nozzle_config)

    init_conditions = InitialConditions(chamber_pressure=PressureValue(1.0, PressureUnit.ATM),
                                        mass_ox=MassValue(6.9, MassUnit.KILOGRAMS),
                                        tank_temp=TemperatureValue(293.1, TemperatureUnit.KELVIN))

    output = sim_loop(sim_config, init_conditions)

    burn_time = output.time[-1]
    chamber_pressure_curve = [t.get_as(PressureUnit.PA) for t in output.chamber_pressure]
    oxidizer_mass_curve = [t.get_as(MassUnit.KILOGRAMS) for t in output.mass_ox]
    fuel_mass_curve = [t.get_as(MassUnit.KILOGRAMS) for t in output.mass_fuel]
    thrust_curve = [t.get_as(MassUnit.NEWTONS) for t in output.thrust]
    
    # plt.plot(output.time, thrust_curve)
    # plt.show()
    print("HRAP Finished")
    return burn_time, chamber_pressure_curve, oxidizer_mass_curve, fuel_mass_curve, thrust_curve

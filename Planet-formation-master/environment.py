"""
环境计算类，用于计算行星的环境参数
"""

import math
from constants import *

class Environment:
    @staticmethod
    def luminosity(mass_ratio):
        """计算恒星光度"""
        if mass_ratio < 1.0:
            n = 1.75 * (mass_ratio - 0.1) + 3.325
        else:
            n = 0.5 * (2.0 - mass_ratio) + 4.4
        return mass_ratio ** n

    @staticmethod
    def orbital_zone(system, orbital_radius):
        """确定行星轨道区域"""
        if orbital_radius < (4.0 * math.sqrt(system.stellar_luminosity_ratio)):
            return 1
        if (orbital_radius >= (4.0 * math.sqrt(system.stellar_luminosity_ratio)) and 
            orbital_radius < (15.0 * math.sqrt(system.stellar_luminosity_ratio))):
            return 2
        return 3

    @staticmethod
    def volume_radius(mass, density):
        """计算行星体积半径"""
        mass = mass * SOLAR_MASS_IN_GRAMS
        volume = mass / density
        return (math.pow((3.0 * volume) / (4.0 * math.pi), (1.0/3.0)) / CM_PER_KM)

    @staticmethod
    def kothari_radius(mass, orbital_radius, giant, zone):
        """使用Kothari公式计算行星半径"""
        if zone == 1:
            if giant:
                atomic_weight = 9.5
                atomic_num = 4.5
            else:
                atomic_weight = 15.0
                atomic_num = 8.0
        elif zone == 2:
            if giant:
                atomic_weight = 2.47
                atomic_num = 2.0
            else:
                atomic_weight = 10.0
                atomic_num = 5.0
        else:
            if giant:
                atomic_weight = 7.0
                atomic_num = 4.0
            else:
                atomic_weight = 10.0
                atomic_num = 5.0

        temp = atomic_weight * atomic_num
        temp = (2.0 * BETA_20 * math.pow(SOLAR_MASS_IN_GRAMS, (1.0/3.0))) / (A1_20 * math.pow(temp, (1.0/3.0)))
        temp2 = A2_20 * math.pow(atomic_weight, (4.0/3.0)) * math.pow(SOLAR_MASS_IN_GRAMS, (2.0/3.0))
        temp2 = temp2 * math.pow(mass, (2.0/3.0))
        temp2 = temp2 / (A1_20 * math.pow(atomic_num, 2.0))
        temp2 = 1.0 + temp2
        temp = temp / temp2
        temp = (temp * math.pow(mass, (1.0/3.0))) / CM_PER_KM
        return temp

    @staticmethod
    def empirical_density(system, mass, orbital_radius, gas_giant):
        """计算行星密度"""
        temp = math.pow(mass * EARTH_MASSES_PER_SOLAR_MASS, (1.0/8.0))
        temp = temp * math.pow(system.r_ecosphere/orbital_radius, (1.0/4.0))
        if gas_giant:
            return temp * 1.2
        else:
            return temp * 5.5

    @staticmethod
    def volume_density(mass, equatorial_radius):
        """计算体积密度"""
        mass = mass * SOLAR_MASS_IN_GRAMS
        equatorial_radius = equatorial_radius * CM_PER_KM
        volume = (4.0 * math.pi * math.pow(equatorial_radius, 3.0)) / 3.0
        return mass / volume

    @staticmethod
    def period(separation, small_mass, large_mass):
        """计算轨道周期"""
        period_in_years = math.sqrt(math.pow(separation, 3.0) / (small_mass + large_mass))
        return period_in_years * DAYS_IN_A_YEAR

    @staticmethod
    def day_length(system, mass, radius, orbital_period, eccentricity, giant):
        """计算日长"""
        system.spin_resonance = False
        k2 = 0.24 if giant else 0.33
        planetary_mass_in_grams = mass * SOLAR_MASS_IN_GRAMS
        equatorial_radius_in_cm = radius * CM_PER_KM
        base_angular_velocity = math.sqrt(2.0 * J * planetary_mass_in_grams / 
                                        (k2 * math.pow(equatorial_radius_in_cm, 2.0)))
        
        change_in_angular_velocity = 0.0
        temp = base_angular_velocity + (change_in_angular_velocity * system.age)
        temp = 1.0 / ((temp / (2.0 * math.pi)) * SECONDS_PER_HOUR)
        
        if temp < orbital_period:
            return temp
            
        spin_resonance_period = ((1.0 - eccentricity) / (1.0 + eccentricity)) * orbital_period
        if eccentricity > 0.01:
            system.callback("...共振...\n")
            temp = spin_resonance_period
            system.spin_resonance = True
        else:
            temp = orbital_period
        return temp

    @staticmethod
    def inclination(system, orbital_radius):
        """计算轨道倾角"""
        temp = math.pow(orbital_radius, 0.2) * system.random_about(EARTH_AXIAL_TILT, 0.4)
        return temp

    @staticmethod
    def escape_velocity(mass, radius):
        """计算逃逸速度"""
        mass_in_grams = mass * SOLAR_MASS_IN_GRAMS
        radius_in_cm = radius * CM_PER_KM
        return math.sqrt(2.0 * GRAV_CONSTANT * mass_in_grams / radius_in_cm)

    @staticmethod
    def rms_velocity(molecular_weight, orbital_radius):
        """计算均方根速度"""
        exospheric_temp = EARTH_EXOSPHERE_TEMP / math.pow(orbital_radius, 2.0)
        return math.sqrt((3.0 * MOLAR_GAS_CONST * exospheric_temp) / molecular_weight) * CM_PER_METER

    @staticmethod
    def molecule_limit(orbital_radius, mass, equatorial_radius):
        """计算可保留的最小分子量"""
        escape_velocity = Environment.escape_velocity(mass, equatorial_radius)
        return ((3.0 * math.pow(GAS_RETENTION_THRESHOLD * CM_PER_METER, 2.0) * 
                MOLAR_GAS_CONST * EARTH_EXOSPHERE_TEMP) / math.pow(escape_velocity, 2.0))

    @staticmethod
    def acceleration(mass, radius):
        """计算表面加速度"""
        return (GRAV_CONSTANT * (mass * SOLAR_MASS_IN_GRAMS) / 
                math.pow(radius * CM_PER_KM, 2.0))

    @staticmethod
    def gravity(acceleration):
        """计算表面重力"""
        return acceleration / EARTH_ACCELERATION

    @staticmethod
    def greenhouse(zone, orbital_radius, greenhouse_radius):
        """计算温室效应"""
        return (orbital_radius < greenhouse_radius) and (zone == 1)

    @staticmethod
    def vol_inventory(system, mass, escape_vel, rms_vel, stellar_mass, zone, greenhouse_effect):
        """计算挥发性气体库存"""
        velocity_ratio = escape_vel / rms_vel
        if velocity_ratio < GAS_RETENTION_THRESHOLD:
            return 0.0
            
        if zone == 1:
            proportion_const = 100000.0
        elif zone == 2:
            proportion_const = 75000.0
        elif zone == 3:
            proportion_const = 250.0
        else:
            proportion_const = 10.0
            system.callback("Error: orbital zone not initialized correctly!\n")
            
        mass_in_earth_units = mass * EARTH_MASSES_PER_SOLAR_MASS
        temp1 = (proportion_const * mass_in_earth_units) / stellar_mass
        temp2 = system.random_about(temp1, 0.2)
        
        if greenhouse_effect:
            return temp2
        return temp2 / 100.0

    @staticmethod
    def pressure(volatile_gas_inventory, equatorial_radius, gravity):
        """计算表面压力"""
        equatorial_radius = EARTH_RADIUS_IN_KM / equatorial_radius
        return volatile_gas_inventory * gravity / math.pow(equatorial_radius, 2.0)

    @staticmethod
    def boiling_point(surface_pressure):
        """计算水的沸点"""
        surface_pressure_in_bars = surface_pressure / MILLIBARS_PER_BAR
        return 1.0 / (math.log(surface_pressure_in_bars) / -5050.5 + 1.0 / 373.0)

    @staticmethod
    def hydrosphere_fraction(volatile_gas_inventory, planetary_radius):
        """计算水圈比例"""
        temp = (0.71 * volatile_gas_inventory / 1000.0) * math.pow(EARTH_RADIUS_IN_KM / planetary_radius, 2.0)
        return min(1.0, temp)

    @staticmethod
    def cloud_fraction(surface_temp, smallest_mw_retained, equatorial_radius, hydrosphere_fraction):
        """计算云层覆盖比例"""
        if smallest_mw_retained > WATER_VAPOR:
            return 0.0
            
        surface_area = 4.0 * math.pi * math.pow(equatorial_radius, 2.0)
        hydrosphere_mass = hydrosphere_fraction * surface_area * EARTH_WATER_MASS_PER_AREA
        water_vapor_in_kg = (0.00000001 * hydrosphere_mass) * math.exp(Q2_36 * (surface_temp - 288.0))
        fraction = CLOUD_COVERAGE_FACTOR * water_vapor_in_kg / surface_area
        return min(1.0, fraction)

    @staticmethod
    def ice_fraction(hydrosphere_fraction, surface_temp):
        """计算冰层覆盖比例"""
        if surface_temp > 328.0:
            surface_temp = 328.0
            
        temp = math.pow(((328.0 - surface_temp) / 90.0), 5.0)
        if temp > (1.5 * hydrosphere_fraction):
            temp = 1.5 * hydrosphere_fraction
        return min(1.0, temp)

    @staticmethod
    def effective_temp(ecosphere_radius, orbital_radius, albedo):
        """计算有效温度"""
        return math.sqrt(ecosphere_radius / orbital_radius) * math.pow((1.0 - albedo) / 0.7, 0.25) * EARTH_EFFECTIVE_TEMP

    @staticmethod
    def greenhouse_rise(optical_depth, effective_temp, surface_pressure):
        """计算温室效应升温"""
        convection_factor = EARTH_CONVECTION_FACTOR * math.pow((surface_pressure / EARTH_SURF_PRES_IN_MILLIBARS), 0.25)
        return (math.pow((1.0 + 0.75 * optical_depth), 0.25) - 1.0) * effective_temp * convection_factor

    @staticmethod
    def planet_albedo(system, water_fraction, cloud_fraction, ice_fraction, surface_pressure):
        """计算行星反照率"""
        rock_fraction = 1.0 - water_fraction - ice_fraction
        components = 0.0
        
        if water_fraction > 0.0:
            components += 1.0
        if ice_fraction > 0.0:
            components += 1.0
        if rock_fraction > 0.0:
            components += 1.0
            
        cloud_adjustment = cloud_fraction / components
        
        if rock_fraction >= cloud_adjustment:
            rock_fraction -= cloud_adjustment
        else:
            rock_fraction = 0.0
            
        if water_fraction > cloud_adjustment:
            water_fraction -= cloud_adjustment
        else:
            water_fraction = 0.0
            
        if ice_fraction > cloud_adjustment:
            ice_fraction -= cloud_adjustment
        else:
            ice_fraction = 0.0
            
        cloud_contribution = cloud_fraction * system.random_about(CLOUD_ALBEDO, 0.2)
        
        if surface_pressure == 0.0:
            rock_contribution = rock_fraction * system.random_about(AIRLESS_ROCKY_ALBEDO, 0.3)
        else:
            rock_contribution = rock_fraction * system.random_about(ROCKY_ALBEDO, 0.1)
            
        water_contribution = water_fraction * system.random_about(WATER_ALBEDO, 0.2)
        
        if surface_pressure == 0.0:
            ice_contribution = ice_fraction * system.random_about(AIRLESS_ICE_ALBEDO, 0.4)
        else:
            ice_contribution = ice_fraction * system.random_about(ICE_ALBEDO, 0.1)
            
        return cloud_contribution + rock_contribution + water_contribution + ice_contribution

    @staticmethod
    def opacity(molecular_weight, surface_pressure):
        """计算不透明度"""
        optical_depth = 0.0
        
        if 0.0 <= molecular_weight < 10.0:
            optical_depth += 3.0
        elif 10.0 <= molecular_weight < 20.0:
            optical_depth += 2.34
        elif 20.0 <= molecular_weight < 30.0:
            optical_depth += 1.0
        elif 30.0 <= molecular_weight < 45.0:
            optical_depth += 0.15
        elif 45.0 <= molecular_weight < 100.0:
            optical_depth += 0.05
            
        if surface_pressure >= (70.0 * EARTH_SURF_PRES_IN_MILLIBARS):
            optical_depth *= 8.333
        elif surface_pressure >= (50.0 * EARTH_SURF_PRES_IN_MILLIBARS):
            optical_depth *= 6.666
        elif surface_pressure >= (30.0 * EARTH_SURF_PRES_IN_MILLIBARS):
            optical_depth *= 3.333
        elif surface_pressure >= (10.0 * EARTH_SURF_PRES_IN_MILLIBARS):
            optical_depth *= 2.0
        elif surface_pressure >= (5.0 * EARTH_SURF_PRES_IN_MILLIBARS):
            optical_depth *= 1.5
            
        return optical_depth

    @staticmethod
    def iterate_surface_temp(system, planet):
        """迭代计算表面温度"""
        albedo = 0.0
        water = 0.0
        clouds = 0.0
        ice = 0.0

        optical_depth = Environment.opacity(planet.molecule_weight, planet.surface_pressure)
        effective_temp = Environment.effective_temp(system.r_ecosphere, planet.a, EARTH_ALBEDO)
        greenhouse_rise = Environment.greenhouse_rise(optical_depth, effective_temp, planet.surface_pressure)
        surface_temp = effective_temp + greenhouse_rise
        previous_temp = surface_temp - 5.0

        while abs(surface_temp - previous_temp) > 1.0:
            previous_temp = surface_temp
            water = Environment.hydrosphere_fraction(planet.volatile_gas_inventory, planet.radius)
            clouds = Environment.cloud_fraction(surface_temp, planet.molecule_weight, planet.radius, water)
            ice = Environment.ice_fraction(water, surface_temp)
            
            if surface_temp >= planet.boil_point or surface_temp <= FREEZING_POINT_OF_WATER:
                water = 0.0
                
            albedo = Environment.planet_albedo(system, water, clouds, ice, planet.surface_pressure)
            optical_depth = Environment.opacity(planet.molecule_weight, planet.surface_pressure)
            effective_temp = Environment.effective_temp(system.r_ecosphere, planet.a, albedo)
            greenhouse_rise = Environment.greenhouse_rise(optical_depth, effective_temp, planet.surface_pressure)
            surface_temp = effective_temp + greenhouse_rise

        planet.hydrosphere = water
        planet.cloud_cover = clouds
        planet.ice_cover = ice
        planet.albedo = albedo
        planet.surface_temp = surface_temp 
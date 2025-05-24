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
        temp = math.pow(orbital_radius, 0.2) * system.random.uniform(EARTH_AXIAL_TILT * 0.6, 
                                                                    EARTH_AXIAL_TILT * 1.4)
        return int(temp) % 360

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
        if orbital_radius < greenhouse_radius:
            return True
        return False

    @staticmethod
    def vol_inventory(system, mass, escape_vel, rms_vel, stellar_mass, zone, greenhouse_effect):
        """计算挥发性气体库存"""
        # 实现挥发性气体库存计算
        pass

    @staticmethod
    def pressure(volatile_gas_inventory, equatorial_radius, gravity):
        """计算表面压力"""
        # 实现表面压力计算
        pass

    @staticmethod
    def boiling_point(surface_pressure):
        """计算水的沸点"""
        # 实现沸点计算
        pass

    @staticmethod
    def hydrosphere_fraction(volatile_gas_inventory, planetary_radius):
        """计算水圈比例"""
        # 实现水圈比例计算
        pass

    @staticmethod
    def cloud_fraction(surface_temp, smallest_mw_retained, equatorial_radius, hydrosphere_fraction):
        """计算云层覆盖比例"""
        # 实现云层覆盖比例计算
        pass

    @staticmethod
    def ice_fraction(hydrosphere_fraction, surface_temp):
        """计算冰层覆盖比例"""
        # 实现冰层覆盖比例计算
        pass

    @staticmethod
    def effective_temp(ecosphere_radius, orbital_radius, albedo):
        """计算有效温度"""
        # 实现有效温度计算
        pass

    @staticmethod
    def greenhouse_rise(optical_depth, effective_temp, surface_pressure):
        """计算温室效应升温"""
        # 实现温室效应升温计算
        pass

    @staticmethod
    def planet_albedo(system, water_fraction, cloud_fraction, ice_fraction, surface_pressure):
        """计算行星反照率"""
        # 实现行星反照率计算
        pass

    @staticmethod
    def opacity(molecular_weight, surface_pressure):
        """计算不透明度"""
        # 实现不透明度计算
        pass

    @staticmethod
    def iterate_surface_temp(system, planet):
        """迭代计算表面温度"""
        # 实现表面温度迭代计算
        pass 
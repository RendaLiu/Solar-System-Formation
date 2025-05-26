"""
行星形成模拟主程序
"""

import os
import argparse
import math
from solar_system import SolarSystem
from accretion import Accretion
from environment import Environment
from constants import *

def get_molecule_name(weight):
    """获取分子量对应的气体名称"""
    if weight < MOLECULAR_HYDROGEN:
        return "H2"
    elif weight < HELIUM:
        return "He"
    elif weight < METHANE:
        return "CH4"
    elif weight < AMMONIA:
        return "NH3"
    elif weight < WATER_VAPOR:
        return "H2O"
    elif weight < NEON:
        return "Ne"
    elif weight < MOLECULAR_NITROGEN:
        return "N2"
    elif weight < CARBON_MONOXIDE:
        return "CO"
    elif weight < NITRIC_OXIDE:
        return "NO"
    elif weight < MOLECULAR_OXYGEN:
        return "O2"
    elif weight < HYDROGEN_SULPHIDE:
        return "H2S"
    elif weight < ARGON:
        return "Ar"
    elif weight < CARBON_DIOXIDE:
        return "CO2"
    elif weight < NITROUS_OXIDE:
        return "N2O"
    elif weight < NITROGEN_DIOXIDE:
        return "NO2"
    elif weight < OZONE:
        return "O3"
    elif weight < SULPHUR_DIOXIDE:
        return "SO2"
    elif weight < SULPHUR_TRIOXIDE:
        return "SO3"
    elif weight < KRYPTON:
        return "Kr"
    elif weight < XENON:
        return "Xe"
    else:
        return ""

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='行星形成模拟程序')
    parser.add_argument('-s', '--seed', type=int, default=None, help='随机数种子')
    parser.add_argument('-v', '--verbose', action='store_true', help='显示详细信息')
    parser.add_argument('-f', '--file', type=str, default='New.System', help='输出文件名')
    parser.add_argument('-c', '--count', type=int, default=float('inf'), help='最大行星数量')
    args = parser.parse_args()

    # 创建太阳系系统，与C#版本保持一致
    system = SolarSystem(verbose=args.verbose, moons=False, callback=print)
    
    # 生成行星系统，使用静态方法
    if args.seed is not None:
        SolarSystem.generate(system, args.seed, args.count)
    else:
        SolarSystem.generate(system, count=args.count)

    # 使用完整路径保存结果到文件
    output_path = os.path.join(os.getcwd(), args.file)
    with open(output_path, 'w', encoding='utf-8') as f:
        # 写入恒星信息
        f.write("                         SYSTEM  CHARACTERISTICS\n")
        f.write(f"Mass of central star:          {system.stellar_mass_ratio:.3f} solar masses\n")
        f.write(f"Luminosity of central star:    {system.stellar_luminosity_ratio:.3f} (relative to the sun)\n")
        f.write(f"Total main sequence lifetime:  {system.main_seq_life/1.0E6:.0f} million years\n")
        f.write(f"Current age of stellar system: {system.age/1.0E6:.0f} million years\n")
        f.write(f"Radius of habitable ecosphere: {system.r_ecosphere:.3f} AU\n\n")

        # 写入行星信息
        for i, planet in enumerate(system.Bodies, 1):  # 使用Bodies而不是bodies
            f.write(f"Planet #{i}:\n")
            
            # 基本信息
            if planet.gas_giant:
                f.write("Gas giant...\n")
            if planet.resonant_period:
                f.write("In resonant period with primary.\n")

            # 轨道参数
            f.write(f"   Distance from primary star (in A.U.): {planet.a:.3f}\n")
            f.write(f"   Eccentricity of orbit:                {planet.e:.3f}\n")
            
            # 物理参数
            f.write(f"   Mass (in Earth masses):               {planet.mass * EARTH_MASSES_PER_SOLAR_MASS:.3f}\n")
            f.write(f"   Equatorial radius (in Km):            {planet.radius:.1f}\n")
            f.write(f"   Density (in g/cc):                    {planet.density:.3f}\n")
            f.write(f"   Escape Velocity (in km/sec):          {planet.escape_velocity/CM_PER_KM:.2f}\n")
            
            # 分子量信息
            molecule_name = get_molecule_name(planet.molecule_weight)
            f.write(f"   Smallest molecular weight retained:   {planet.molecule_weight:.2f}")
            if molecule_name:
                f.write(f"   ({molecule_name})\n")
            else:
                f.write("\n")

            # 表面条件
            f.write(f"   Surface acceleration (in cm/sec2):    {planet.surface_accel:.2f}\n")
            if not planet.gas_giant:
                f.write(f"   Surface Gravity (in Earth gees):      {planet.surface_grav:.2f}\n")
                if planet.boil_point > 0.1:
                    f.write(f"   Boiling point of water (celcius):     {planet.boil_point - KELVIN_CELCIUS_DIFFERENCE:.1f}\n")
                if planet.surface_pressure > 0.00001:
                    f.write(f"   Surface Pressure (in atmospheres):    {planet.surface_pressure/1000.0:.3f}")
                    f.write("     RUNAWAY GREENHOUSE EFFECT\n" if planet.greenhouse_effect else "\n")
                f.write(f"   Surface temperature (Celcius):        {planet.surface_temp - KELVIN_CELCIUS_DIFFERENCE:.2f}\n")
                
                # 水、云、冰覆盖率
                if planet.hydrosphere > 0.01:
                    f.write(f"   Hydrosphere percentage:               {planet.hydrosphere * 100:.2f}\n")
                if planet.cloud_cover > 0.01:
                    f.write(f"   Cloud cover percentage:               {planet.cloud_cover * 100:.2f}\n")
                if planet.ice_cover > 0.01:
                    f.write(f"   Ice cover percentage:                 {planet.ice_cover * 100:.2f}\n")

            # 其他参数
            f.write(f"   Axial tilt (in degrees):              {planet.axial_tilt}\n")
            f.write(f"   Planetary albedo:                     {planet.albedo:.3f}\n")
            f.write(f"   Length of year (in years):            {planet.orbital_period/365.25:.2f}\n")
            f.write(f"   Length of day (in hours):             {planet.day:.2f}\n\n")

    print()
    print(f'Done! System definition written to "{output_path}"')

if __name__ == "__main__":
    main() 
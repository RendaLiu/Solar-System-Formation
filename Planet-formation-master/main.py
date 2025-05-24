"""
行星形成模拟主程序
"""

from solar_system import SolarSystem
from accretion import Accretion
from environment import Environment
import argparse
import math
from constants import *

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='行星形成模拟程序')
    parser.add_argument('-s', '--seed', type=int, default=None, help='随机数种子')
    parser.add_argument('-v', '--verbose', action='store_true', help='显示详细信息')
    parser.add_argument('-f', '--file', type=str, default='New.System', help='输出文件名')
    parser.add_argument('-c', '--count', type=int, default=None, help='最大行星数量')
    args = parser.parse_args()

    # 创建太阳系系统
    system = SolarSystem(verbose=args.verbose, use_seed=args.seed is not None)
    if args.seed is not None:
        system.random.seed(args.seed)

    # 设置恒星参数
    system.set_stellar_parameters(1.0, 1.0)  # 太阳质量比和光度比

    # 计算尘埃云边界
    inner_dust = Accretion.innermost_planet(system.stellar_mass_ratio)
    outer_dust = Accretion.outermost_planet(system.stellar_mass_ratio)

    # 初始化尘埃云
    Accretion.set_initial_conditions(system, inner_dust, outer_dust)

    # 开始行星形成过程
    while system.dust_left:
        # 随机选择位置注入初始质量
        a = system.random.uniform(inner_dust, outer_dust)
        e = system.random.uniform(0.0, 0.1)  # 初始离心率
        mass = PROTOPLANET_MASS

        if system.verbose:
            system.callback(f"检查 {a:.3f} AU 位置\n")

        # 检查是否有可用尘埃
        if not Accretion.dust_available(system, 
                                      Accretion.inner_effect_limit(system, a, e, mass),
                                      Accretion.outer_effect_limit(system, a, e, mass)):
            if system.verbose:
                system.callback(".. 失败\n")
            continue

        if system.verbose:
            system.callback(".. 注入原行星\n")

        # 计算尘埃密度
        system.dust_density = (DUST_DENSITY_COEFF * 
                             math.sqrt(system.stellar_mass_ratio) * 
                             math.exp(-ALPHA * math.pow(a, 1.0/N)))

        # 计算临界质量
        crit_mass = Accretion.critical_limit(a, e, system.stellar_luminosity_ratio)

        # 执行吸积过程
        mass = Accretion.accrete_dust(system, mass, a, e, crit_mass, 
                                    inner_dust, outer_dust)

        # 处理星子合并
        if mass != 0.0 and mass != PROTOPLANET_MASS:
            Accretion.coalesce_planetesimals(system, a, e, mass, crit_mass,
                                           system.stellar_luminosity_ratio,
                                           inner_dust, outer_dust)
        else:
            system.callback(".. 由于大质量邻居而失败\n")

        # 检查是否达到最大行星数量
        if args.count is not None and len(system.bodies) >= args.count:
            break

    # 计算行星环境参数
    for planet in system.bodies:
        # 计算基本物理参数
        zone = Environment.orbital_zone(system, planet.a)
        planet.radius = Environment.kothari_radius(planet.mass, planet.a, 
                                                 planet.gas_giant, zone)
        planet.density = Environment.empirical_density(system, planet.mass, 
                                                     planet.a, planet.gas_giant)
        
        # 计算环境参数
        planet.escape_velocity = Environment.escape_velocity(planet.mass, planet.radius)
        planet.molecule_weight = Environment.molecule_limit(planet.a, planet.mass, 
                                                          planet.radius)
        planet.surface_accel = Environment.acceleration(planet.mass, planet.radius)
        planet.surface_grav = Environment.gravity(planet.surface_accel)
        
        # 计算轨道参数
        planet.orbital_period = Environment.period(planet.a, planet.mass, 
                                                 system.stellar_mass_ratio)
        planet.day = Environment.day_length(system, planet.mass, planet.radius,
                                          planet.orbital_period, planet.e, 
                                          planet.gas_giant)
        planet.axial_tilt = Environment.inclination(system, planet.a)

    # 输出系统信息
    system.print_system_info()

    # 保存结果到文件
    with open(args.file, 'w') as f:
        f.write("                         SYSTEM  CHARACTERISTICS\n")
        f.write(f"Mass of central star:          {system.stellar_mass_ratio:.3f} solar masses\n")
        f.write(f"Luminosity of central star:    {system.stellar_luminosity_ratio:.3f} (relative to the sun)\n")
        f.write(f"Total main sequence lifetime:  {system.main_seq_life/1.0E6:.0f} million years\n")
        f.write(f"Current age of stellar system: {system.age/1.0E6:.0f} million years\n")
        f.write(f"Radius of habitable ecosphere: {system.r_ecosphere:.3f} AU\n\n")

        for i, planet in enumerate(system.bodies, 1):
            f.write(f"Planet #{i}:\n")
            if planet.gas_giant:
                f.write("Gas giant...\n")
            f.write(f"   Distance from primary star (in A.U.): {planet.a:.3f}\n")
            f.write(f"   Eccentricity of orbit:                {planet.e:.3f}\n")
            f.write(f"   Mass (in Earth masses):               {planet.mass * 332775.64:.3f}\n")
            f.write(f"   Equatorial radius (in Km):            {planet.radius:.1f}\n")
            f.write(f"   Density (in g/cc):                    {planet.density:.3f}\n")
            f.write(f"   Escape Velocity (in km/sec):          {planet.escape_velocity:.2f}\n")
            f.write(f"   Surface acceleration (in cm/sec2):    {planet.surface_accel:.2f}\n")
            f.write(f"   Surface Gravity (in Earth gees):      {planet.surface_grav:.2f}\n")
            f.write(f"   Length of year (in years):            {planet.orbital_period/365.25:.2f}\n")
            f.write(f"   Length of day (in hours):             {planet.day:.2f}\n\n")

if __name__ == "__main__":
    main() 
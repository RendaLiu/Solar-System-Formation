"""
太阳系类，用于管理整个行星系统
"""

import random
import math
from dust import Dust
from planet import Planet
from stellartype import StellarType
from environment import Environment
from constants import *
from accretion import Accretion

class SolarSystem:
    def __init__(self, verbose=False, moons=False, callback=None):
        # 恒星参数
        self.radians_per_rotation = 2.0 * math.pi
        self.stellar_mass_ratio = 1.0
        self.stellar_luminosity_ratio = 1.0
        self.stellar_radius_ratio = 1.0
        self.stellar_temp = 0.0
        self.main_seq_life = 1.0e10
        self.age = 4.6e9
        self.r_ecosphere = 1.0
        self.r_greenhouse = 1.0
        self.type = None
        
        # 行星系统参数
        self.planet_head = None
        self.bodies = []  # 改为普通列表
        self.Bodies = None  # 添加数组形式属性
        self.anum = 0.0
        self.spin_resonance = False
        self.moons = moons
        
        # 尘埃云参数
        self.dust_left = True
        self.dust_head = None
        self.dust_density = 0.0
        self.cloud_eccentricity = 0.2
        
        # 计算参数
        self.r_inner = 0.0
        self.r_outer = 0.0
        self.reduced_mass = 0.0
        
        # 系统参数
        self.verbose = verbose
        self.random = random.Random()
        self.callback = callback if callback else print

    @staticmethod
    def generate(system, seed=None, count=float('inf')):
        """静态生成方法"""
        if seed is not None:
            system.random = random.Random(seed)
        system._generate_internal(count)

    def random_about(self, value, variation):
        """在指定值附近生成随机数"""
        return value * (1.0 + (self.random.random() - 0.5) * variation)

    def random_eccentricity(self):
        """生成随机离心率"""
        return self.random.random() * 0.1

    def random_range(self, min_val, max_val):
        """生成指定范围内的随机数"""
        return min_val + (max_val - min_val) * self.random.random()

    def set_stellar_parameters(self, mass_ratio=None, luminosity_ratio=None):
        """设置恒星参数"""
        if mass_ratio is not None and luminosity_ratio is not None:
            self.stellar_mass_ratio = mass_ratio
            self.stellar_luminosity_ratio = luminosity_ratio
        else:
            # 随机生成恒星参数
            self.stellar_mass_ratio = self.random_range(0.6, 1.3)
            self.stellar_radius_ratio = math.floor(
                self.random_about(
                    math.pow(self.stellar_mass_ratio, 1.0/3.0), 
                    0.05
                ) * 1000.0
            ) / 1000.0
            self.stellar_luminosity_ratio = Environment.luminosity(self.stellar_mass_ratio)
            
        # 计算其他参数
        self.stellar_temp = math.floor(
            5650 * math.sqrt(
                math.sqrt(self.stellar_luminosity_ratio) / self.stellar_radius_ratio
            )
        )
        self.main_seq_life = 1.0E10 * (
            self.stellar_mass_ratio / self.stellar_luminosity_ratio
        )
        self.age = self.random_range(
            1.0E9, 
            min(6.0E9, self.main_seq_life)
        )
        self.r_ecosphere = math.sqrt(self.stellar_luminosity_ratio)
        self.r_greenhouse = self.r_ecosphere * GREENHOUSE_EFFECT_CONST
        self.type = StellarType.get_stellar_type_temp(self.stellar_temp)

    def calculate_ecosphere_radius(self):
        """计算宜居带半径"""
        return 1.0 * (self.stellar_luminosity_ratio ** 0.5)
        
    def add_planet(self, planet):
        """添加行星到系统"""
        # 添加到列表并保持有序
        self.bodies.append(planet)
        self.bodies.sort(key=lambda x: x.a)
        
        # 维护链表结构
        if self.planet_head is None:
            self.planet_head = planet
            planet.next_planet = None
        else:
            current = self.planet_head
            if planet.a < current.a:
                planet.next_planet = current
                self.planet_head = planet
            else:
                while current.next_planet and current.next_planet.a < planet.a:
                    current = current.next_planet
                planet.next_planet = current.next_planet
                current.next_planet = planet

    def print_system_info(self):
        """打印系统信息"""
        self.callback(f"恒星质量比: {self.stellar_mass_ratio:.3f}")
        self.callback(f"恒星光度比: {self.stellar_luminosity_ratio:.3f}")
        self.callback(f"系统年龄: {self.age/1e6:.0f} 百万年")
        self.callback(f"宜居带半径: {self.r_ecosphere:.3f} AU\n")
        
        planet = self.planet_head
        counter = 1
        while planet:
            self.callback(f"行星 #{counter}:")
            if planet.gas_giant:
                self.callback("气态巨行星")
            self.callback(f"  距恒星距离 (AU): {planet.a:.3f}")
            self.callback(f"  轨道离心率: {planet.e:.3f}")
            self.callback(f"  质量 (地球质量): {planet.mass * 332775.64:.3f}")
            self.callback(f"  赤道半径 (km): {planet.radius:.1f}")
            self.callback(f"  密度 (g/cm³): {planet.density:.3f}")
            self.callback(f"  逃逸速度 (km/s): {planet.escape_velocity:.2f}")
            self.callback(f"  表面温度 (℃): {planet.surface_temp:.2f}\n")
            counter += 1
            planet = planet.next_planet 

    def _generate_internal(self, count):
        """内部生成实现"""
        # 设置恒星参数
        self.set_stellar_parameters()
        
        # 计算尘埃云边界
        inner_dust = Accretion.innermost_planet(self.stellar_mass_ratio)
        outer_dust = Accretion.stellar_dust_limit(self, self.stellar_mass_ratio)
        
        # 初始化尘埃云
        Accretion.set_initial_conditions(self, inner_dust, outer_dust)
        
        # 生成行星
        i = 0
        planet = Accretion.distribute_planetary_masses(
            self, 
            self.stellar_mass_ratio,
            self.stellar_luminosity_ratio,
            inner_dust,
            outer_dust
        )
        
        while planet and i < count:
            # 计算轨道带
            planet.orbit_zone = Environment.orbital_zone(self, planet.a)
            
            # 计算物理参数
            if planet.gas_giant:
                planet.density = Environment.empirical_density(
                    self, planet.mass, planet.a, planet.gas_giant
                )
                planet.radius = Environment.volume_radius(planet.mass, planet.density)
            else:
                planet.radius = Environment.kothari_radius(
                    planet.mass, planet.a, planet.gas_giant, planet.orbit_zone
                )
                planet.density = Environment.volume_density(planet.mass, planet.radius)
            
            # 计算轨道参数
            self._calculate_orbital_parameters(planet)
            
            # 计算大气参数
            self._calculate_atmospheric_parameters(planet)
            
            # 生成卫星系统
            if self.moons:
                self._generate_moons(planet)
            
            # 继续处理下一个行星
            if i + 1 < count:
                planet = planet.next_planet
                i += 1
            else:
                planet.next_planet = None
                # 移除超出范围的行星
                self.bodies = [p for p in self.bodies if p.a <= planet.a]
                self.bodies.sort(key=lambda x: x.a)
                break
        
        # 更新Bodies数组
        self.Bodies = list(self.bodies)

    def _calculate_orbital_parameters(self, planet):
        """计算轨道参数"""
        planet.orbital_period = Environment.period(
            planet.a, planet.mass, self.stellar_mass_ratio
        )
        planet.day = Environment.day_length(
            self, planet.mass, planet.radius,
            planet.orbital_period, planet.e, planet.gas_giant
        )
        planet.resonant_period = self.spin_resonance
        planet.axial_tilt = Environment.inclination(self, planet.a)

    def _calculate_atmospheric_parameters(self, planet):
        """计算大气参数"""
        planet.escape_velocity = Environment.escape_velocity(planet.mass, planet.radius)
        planet.surface_accel = Environment.acceleration(planet.mass, planet.radius)
        planet.rms_velocity = Environment.rms_velocity(MOLECULAR_NITROGEN, planet.a)
        planet.molecule_weight = Environment.molecule_limit(
            planet.a, planet.mass, planet.radius
        )
        
        if planet.gas_giant:
            planet.greenhouse_effect = False
            planet.hydrosphere = INCREDIBLY_LARGE_NUMBER
            planet.albedo = self.random_about(GAS_GIANT_ALBEDO, 0.1)
        else:
            self._calculate_surface_conditions(planet)

    def _calculate_surface_conditions(self, planet):
        """计算表面条件"""
        planet.surface_grav = Environment.gravity(planet.surface_accel)
        planet.greenhouse_effect = Environment.greenhouse(
            planet.orbit_zone, planet.a, self.r_greenhouse
        )
        planet.volatile_gas_inventory = Environment.vol_inventory(
            self, planet.mass, planet.escape_velocity,
            planet.rms_velocity, self.stellar_mass_ratio,
            planet.orbit_zone, planet.greenhouse_effect
        )
        planet.surface_pressure = Environment.pressure(
            planet.volatile_gas_inventory,
            planet.radius,
            planet.surface_grav
        )
        planet.boil_point = (
            0.0 if planet.surface_pressure == 0.0
            else Environment.boiling_point(planet.surface_pressure)
        )
        Environment.iterate_surface_temp(self, planet)

    def _generate_moons(self, planet):
        """生成卫星系统"""
        planet.first_moon = Accretion.distribute_moon_masses(
            self, planet.mass, planet.radius
        )
        
        moon = planet.first_moon
        while moon:
            self._calculate_moon_parameters(moon, planet)
            planet.bodies_orbiting.append(moon)
            moon = moon.next_planet
            
        planet.bodies_orbiting = list(planet.bodies_orbiting)

    def _calculate_moon_parameters(self, moon, planet):
        """计算卫星参数"""
        moon.radius = Environment.kothari_radius(
            moon.mass, 0, False, planet.orbit_zone
        )
        moon.density = Environment.volume_density(moon.mass, moon.radius)
        moon.density = max(
            1.5,
            self.random_range(1.5, moon.density * 1.1)
        )
        moon.radius = Environment.volume_radius(moon.mass, moon.density)
        
        # 计算轨道参数
        moon.orbital_period = Environment.period(
            moon.a, moon.mass, planet.mass
        )
        moon.day = Environment.day_length(
            self, moon.mass, moon.radius,
            moon.orbital_period, moon.e, False
        )
        moon.resonant_period = self.spin_resonance
        moon.axial_tilt = Environment.inclination(self, moon.a)
        
        # 计算表面参数
        self._calculate_atmospheric_parameters(moon)
        Environment.iterate_surface_temp(self, planet, moon) 
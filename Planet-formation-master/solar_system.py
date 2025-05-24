"""
太阳系类，用于管理整个行星系统
"""

import random
from dust import Dust
from planet import Planet

class SolarSystem:
    def __init__(self, verbose=False, use_seed=False, callback=None):
        # 恒星参数
        self.stellar_mass_ratio = 1.0  # 恒星质量比（相对于太阳）
        self.stellar_luminosity_ratio = 1.0  # 恒星光度比（相对于太阳）
        self.main_seq_life = 1.0e10  # 主序星寿命（年）
        self.age = 4.6e9  # 系统年龄（年）
        self.r_ecosphere = 1.0  # 宜居带半径（AU）
        
        # 尘埃云参数
        self.dust_head = None  # 第一个尘埃带
        self.dust_left = True  # 是否还有尘埃
        self.dust_density = 0.0  # 尘埃密度
        self.cloud_eccentricity = 0.2  # 尘埃云离心率
        
        # 行星系统参数
        self.planet_head = None  # 第一个行星
        self.bodies = []  # 所有天体的列表
        
        # 计算参数
        self.r_inner = 0.0  # 内影响范围
        self.r_outer = 0.0  # 外影响范围
        self.reduced_mass = 0.0  # 约化质量
        
        # 系统参数
        self.verbose = verbose  # 是否输出详细信息
        self.random = random.Random()  # 随机数生成器
        if use_seed:
            self.random.seed(42)  # 使用固定种子
        self.callback = callback if callback else print  # 回调函数
        
    def set_stellar_parameters(self, mass_ratio, luminosity_ratio):
        """设置恒星参数"""
        self.stellar_mass_ratio = mass_ratio
        self.stellar_luminosity_ratio = luminosity_ratio
        self.r_ecosphere = self.calculate_ecosphere_radius()
        
    def calculate_ecosphere_radius(self):
        """计算宜居带半径"""
        return 1.0 * (self.stellar_luminosity_ratio ** 0.5)
        
    def add_planet(self, planet):
        """添加行星到系统"""
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
        self.bodies.append(planet)
        
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
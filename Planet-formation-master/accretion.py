"""
吸积模拟类，实现行星形成的主要物理过程
"""

import math
from constants import *
from dust import Dust
from planet import Planet

class Accretion:
    @staticmethod
    def set_initial_conditions(system, inner_limit_of_dust, outer_limit_of_dust):
        """设置初始条件"""
        system.dust_head = Dust()
        system.planet_head = None
        system.dust_head.next_band = None
        system.dust_head.outer_edge = outer_limit_of_dust
        system.dust_head.inner_edge = inner_limit_of_dust
        system.dust_head.dust_present = True
        system.dust_head.gas_present = True
        system.dust_left = True
        system.cloud_eccentricity = 0.2

    @staticmethod
    def stellar_dust_limit(system, stellar_mass_ratio):
        """计算恒星的尘埃限制"""
        return 200.0 * (stellar_mass_ratio ** (1.0/3.0))

    @staticmethod
    def innermost_planet(stellar_mass_ratio):
        """计算最内层行星位置"""
        return 0.3 * (stellar_mass_ratio ** (1.0/3.0))

    @staticmethod
    def outermost_planet(stellar_mass_ratio):
        """计算最外层行星位置"""
        return 50.0 * (stellar_mass_ratio ** (1.0/3.0))

    @staticmethod
    def inner_effect_limit(system, a, e, mass):
        """计算内影响范围"""
        return a * (1.0 - e) * (1.0 - mass) / (1.0 + system.cloud_eccentricity)

    @staticmethod
    def outer_effect_limit(system, a, e, mass):
        """计算外影响范围"""
        return a * (1.0 + e) * (1.0 + system.reduced_mass) / (1.0 - system.cloud_eccentricity)

    @staticmethod
    def dust_available(system, inside_range, outside_range):
        """检查指定范围内是否有可用尘埃"""
        current_dust_band = system.dust_head
        while current_dust_band and current_dust_band.outer_edge < inside_range:
            current_dust_band = current_dust_band.next_band
        dust_here = current_dust_band is not None and current_dust_band.dust_present
        while current_dust_band and current_dust_band.inner_edge < outside_range:
            dust_here = dust_here or current_dust_band.dust_present
            current_dust_band = current_dust_band.next_band
        return dust_here

    @staticmethod
    def update_dust_lanes(system, min_range, max_range, mass, crit_mass, body_inner_bound, body_outer_bound):
        """更新尘埃带状态"""
        system.dust_left = False
        gas = not (mass > crit_mass)
        node1 = system.dust_head
        
        while node1:
            if node1.inner_edge < min_range and node1.outer_edge > max_range:
                # 创建新的尘埃带
                node2 = Dust()
                node2.inner_edge = min_range
                node2.outer_edge = max_range
                node2.gas_present = gas if node1.gas_present else False
                node2.dust_present = False
                
                node3 = Dust()
                node3.inner_edge = max_range
                node3.outer_edge = node1.outer_edge
                node3.gas_present = node1.gas_present
                node3.dust_present = node1.dust_present
                node3.next_band = node1.next_band
                
                node1.next_band = node2
                node2.next_band = node3
                node1.outer_edge = min_range
                node1 = node3.next_band
            else:
                node1 = node1.next_band

    @staticmethod
    def collect_dust(system, last_mass, a, e, crit_mass, dust_band):
        """收集尘埃物质"""
        if dust_band is None:
            return 0.0
            
        temp = last_mass / (1.0 + last_mass)
        system.reduced_mass = temp ** (1.0/4.0)
        system.r_inner = Accretion.inner_effect_limit(system, a, e, system.reduced_mass)
        system.r_outer = Accretion.outer_effect_limit(system, a, e, system.reduced_mass)
        
        if system.r_inner < 0.0:
            system.r_inner = 0.0
            
        temp_density = 0.0 if not dust_band.dust_present else system.dust_density
        if last_mass < crit_mass or not dust_band.gas_present:
            mass_density = temp_density
        else:
            mass_density = K * temp_density / (1.0 + math.sqrt(crit_mass/last_mass) * (K - 1.0))
            
        if dust_band.outer_edge <= system.r_inner or dust_band.inner_edge >= system.r_outer:
            return Accretion.collect_dust(system, last_mass, a, e, crit_mass, dust_band.next_band)
            
        bandwidth = system.r_outer - system.r_inner
        temp1 = max(0.0, system.r_outer - dust_band.outer_edge)
        width = bandwidth - temp1
        temp2 = max(0.0, dust_band.inner_edge - system.r_inner)
        width = width - temp2
        
        temp = 4.0 * math.pi * (a**2) * system.reduced_mass * (1.0 - e * (temp1 - temp2) / bandwidth)
        volume = temp * width
        
        return volume * mass_density + Accretion.collect_dust(system, last_mass, a, e, crit_mass, dust_band.next_band)

    @staticmethod
    def critical_limit(orbital_radius, eccentricity, stellar_luminosity_ratio):
        """计算临界质量"""
        perihelion_dist = orbital_radius - orbital_radius * eccentricity
        temp = perihelion_dist * math.sqrt(stellar_luminosity_ratio)
        return B * (temp ** -0.75)

    @staticmethod
    def accrete_dust(system, seed_mass, a, e, crit_mass, body_inner_bound, body_outer_bound):
        """执行尘埃吸积过程"""
        temp_mass = 0.0
        new_mass = seed_mass
        
        while True:
            temp_mass = new_mass
            new_mass = Accretion.collect_dust(system, new_mass, a, e, crit_mass, system.dust_head)
            if new_mass - temp_mass < 0.0001 * temp_mass:
                break
                
        seed_mass += new_mass
        Accretion.update_dust_lanes(system, system.r_inner, system.r_outer, 
                                  seed_mass, crit_mass, body_inner_bound, body_outer_bound)
        return seed_mass

    @staticmethod
    def coalesce_planetesimals(system, a, e, mass, crit_mass, stellar_luminosity_ratio, 
                              body_inner_bound, body_outer_bound):
        """模拟星子合并过程"""
        coalesced = False
        node1 = system.planet_head
        
        while node1:
            temp = node1.a - a
            if temp > 0.0:
                dist1 = a * (1.0 + e) * (1.0 + system.reduced_mass) - a
                system.reduced_mass = (node1.mass / (1.0 + node1.mass)) ** (1.0/4.0)
                dist2 = node1.a - node1.a * (1.0 - node1.e) * (1.0 - system.reduced_mass)
            else:
                dist1 = a - a * (1.0 - e) * (1.0 - system.reduced_mass)
                system.reduced_mass = (node1.mass / (1.0 + node1.mass)) ** (1.0/4.0)
                dist2 = node1.a * (1.0 + node1.e) * (1.0 + system.reduced_mass) - node1.a
                
            if abs(temp) <= abs(dist1) or abs(temp) <= abs(dist2):
                system.callback("星子碰撞！\n")
                a3 = (node1.mass + mass) / (node1.mass/node1.a + mass/a)
                temp = node1.mass * math.sqrt(node1.a) * math.sqrt(1.0 - node1.e**2)
                temp += mass * math.sqrt(a) * math.sqrt(1.0 - e**2)
                temp = temp / ((node1.mass + mass) * math.sqrt(a3))
                temp = 1.0 - temp**2
                if temp < 0.0 or temp >= 1.0:
                    temp = 0.0
                e = math.sqrt(temp)
                temp = node1.mass + mass
                temp = Accretion.accrete_dust(system, temp, a3, e, stellar_luminosity_ratio,
                                           body_inner_bound, body_outer_bound)
                node1.a = a3
                node1.e = e
                node1.mass = temp
                coalesced = True
                break
            node1 = node1.next_planet
            
        if not coalesced:
            new_planet = Planet()
            new_planet.a = a
            new_planet.e = e
            new_planet.mass = mass
            new_planet.gas_giant = mass >= crit_mass
            system.add_planet(new_planet) 
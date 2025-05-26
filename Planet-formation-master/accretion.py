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

    @staticmethod #good
    def stellar_dust_limit(system, stellar_mass_ratio):
        """计算恒星的尘埃限制"""
        return 200.0 * (stellar_mass_ratio ** (1.0/3.0))

    @staticmethod #good
    def innermost_planet(stellar_mass_ratio):
        """计算最内层行星位置"""
        return 0.3 * (stellar_mass_ratio ** (1.0/3.0))

    @staticmethod #good
    def outermost_planet(stellar_mass_ratio):
        """计算最外层行星位置"""
        return 50.0 * (stellar_mass_ratio ** (1.0/3.0))

    @staticmethod #good
    def inner_effect_limit(system, a, e, mass):
        """计算内影响范围"""
        return a * (1.0 - e) * (1.0 - mass) / (1.0 + system.cloud_eccentricity)

    @staticmethod #good
    def outer_effect_limit(system, a, e, mass):
        """计算外影响范围"""
        return a * (1.0 + e) * (1.0 + system.reduced_mass) / (1.0 - system.cloud_eccentricity)

    @staticmethod #good
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

    @staticmethod #good
    def update_dust_lanes(system, min_range, max_range, mass, crit_mass, body_inner_bound, body_outer_bound):
        """更新尘埃带状态"""
        # 初始化状态
        system.dust_left = False
        gas = not (mass > crit_mass)
        node1 = system.dust_head
        
        # 第一轮：更新尘埃带
        while node1:
            if node1.inner_edge < min_range and node1.outer_edge > max_range:
                # 情况1：尘埃带完全包含在范围内
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
                
            elif node1.inner_edge < max_range and node1.outer_edge > max_range:
                # 情况2：尘埃带与范围右边界相交
                node2 = Dust()
                node2.next_band = node1.next_band
                node2.dust_present = node1.dust_present
                node2.gas_present = node1.gas_present
                node2.outer_edge = node1.outer_edge
                node2.inner_edge = max_range
                
                node1.next_band = node2
                node1.outer_edge = max_range
                if node1.gas_present:
                    node1.gas_present = gas
                node1.dust_present = False
                node1 = node2.next_band
                
            elif node1.inner_edge < min_range and node1.outer_edge > min_range:
                # 情况3：尘埃带与范围左边界相交
                node2 = Dust()
                node2.next_band = node1.next_band
                node2.dust_present = False
                if node1.gas_present:
                    node2.gas_present = gas
                node2.outer_edge = node1.outer_edge
                node2.inner_edge = min_range
                
                node1.next_band = node2
                node1.outer_edge = min_range
                node1 = node2.next_band
                
            elif node1.inner_edge >= min_range and node1.outer_edge <= max_range:
                # 情况4：尘埃带完全在范围内
                if node1.gas_present:
                    node1.gas_present = gas
                node1.dust_present = False
                node1 = node1.next_band
                
            elif node1.outer_edge < min_range or node1.inner_edge > max_range:
                # 情况5：尘埃带与范围无交集
                node1 = node1.next_band
        
        # 第二轮：检查剩余尘埃并合并相邻的相同状态尘埃带
        node1 = system.dust_head
        while node1:
            # 检查是否有尘埃剩余
            if (node1.dust_present and 
                node1.outer_edge >= body_inner_bound and 
                node1.inner_edge <= body_outer_bound):
                system.dust_left = True
                
            node2 = node1.next_band
            # 合并相邻的相同状态尘埃带
            if (node2 and 
                node1.dust_present == node2.dust_present and 
                node1.gas_present == node2.gas_present):
                node1.outer_edge = node2.outer_edge
                node1.next_band = node2.next_band
                node2 = None
                
            node1 = node1.next_band

    @staticmethod #good 递归调用了collect_dust
    def collect_dust(system, last_mass, a, e, crit_mass, dust_band):
        """收集尘埃物质"""
        if dust_band is None:
            return 0.0
        
        mass_density = 0.0
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

    @staticmethod #good
    def critical_limit(orbital_radius, eccentricity, stellar_luminosity_ratio):
        """计算临界质量"""
        perihelion_dist = orbital_radius - orbital_radius * eccentricity
        temp = perihelion_dist * math.sqrt(stellar_luminosity_ratio)
        return B * (temp ** -0.75)

    @staticmethod #good
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
                # 远日点距离
                dist1 = a * (1.0 + e) * (1.0 + system.reduced_mass) - a
                system.reduced_mass = (node1.mass / (1.0 + node1.mass)) ** (1.0/4.0)
                dist2 = node1.a - node1.a * (1.0 - node1.e) * (1.0 - system.reduced_mass)
            else:
                # 近日点距离
                dist1 = a - a * (1.0 - e) * (1.0 - system.reduced_mass)
                system.reduced_mass = (node1.mass / (1.0 + node1.mass)) ** (1.0/4.0)
                dist2 = node1.a * (1.0 + node1.e) * (1.0 + system.reduced_mass) - node1.a
                
            if abs(temp) <= abs(dist1) or abs(temp) <= abs(dist2):
                system.callback("星子碰撞！\n")
                # 计算新的轨道半径
                a3 = (node1.mass + mass) / (node1.mass/node1.a + mass/a)
                # 计算新的偏心率
                temp = node1.mass * math.sqrt(node1.a) * math.sqrt(1.0 - node1.e**2)
                # 修正：添加额外的sqrt
                temp += mass * math.sqrt(a) * math.sqrt(math.sqrt(1.0 - e**2))
                temp = temp / ((node1.mass + mass) * math.sqrt(a3))
                temp = 1.0 - temp**2
                
                if temp < 0.0 or temp >= 1.0:
                    temp = 0.0
                e = math.sqrt(temp)
                
                # 计算新的质量
                temp = node1.mass + mass
                temp = Accretion.accrete_dust(system, temp, a3, e, stellar_luminosity_ratio,
                                        body_inner_bound, body_outer_bound)
                
                # 更新行星参数
                node1.a = a3
                node1.e = e
                node1.mass = temp
                coalesced = True
                break
                
            node1 = node1.next_planet
        
        if not coalesced:
            # 创建新行星
            new_planet = Planet()
            new_planet.a = a
            new_planet.e = e
            new_planet.mass = mass
            new_planet.gas_giant = mass >= crit_mass
            # 使用solar_system中的add_planet方法添加行星
            system.add_planet(new_planet)

    @staticmethod
    def distribute_planetary_masses(system, stellar_mass_ratio, stellar_luminosity_ratio, inner_dust, outer_dust):
        """分配行星质量"""
        # 设置初始参数
        planetesimal_inner_bound = inner_dust
        planetesimal_outer_bound = outer_dust
        
        while Accretion.dust_available(system, inner_dust, outer_dust):
            # 随机生成位置和质量
            a = system.random_range(planetesimal_inner_bound, planetesimal_outer_bound)
            e = system.random_eccentricity()
            
            if system.verbose:
                system.callback(f"Checking {a} AU.\n")
                
            # 检查范围内是否有可用尘埃
            if not Accretion.dust_available(system, 
                                          Accretion.inner_effect_limit(system, a, e, PROTOPLANET_MASS),
                                          Accretion.outer_effect_limit(system, a, e, PROTOPLANET_MASS)):
                if system.verbose:
                    system.callback(".. failed.\n")
                continue
                
            system.callback(".. Injecting protoplanet.\n")
            
            # 计算尘埃密度和临界质量
            system.dust_density = (DUST_DENSITY_COEFF * 
                                 math.sqrt(stellar_mass_ratio) * 
                                 math.exp(-ALPHA * math.pow(a, 1.0/N)))
            crit_mass = Accretion.critical_limit(a, e, stellar_luminosity_ratio)
            
            # 执行吸积过程
            mass = Accretion.accrete_dust(system, PROTOPLANET_MASS, a, e, crit_mass, 
                                        planetesimal_inner_bound, planetesimal_outer_bound)
            
            # 检查是否需要合并
            if mass != 0.0 and mass != PROTOPLANET_MASS:
                Accretion.coalesce_planetesimals(system, a, e, mass, crit_mass, 
                                               stellar_luminosity_ratio,
                                               planetesimal_inner_bound, 
                                               planetesimal_outer_bound)
            else:
                system.callback(".. failed due to large neighbor.\n")
                
        return system.planet_head

    @staticmethod
    def distribute_moon_masses(system, planetary_mass, plan_radius):
        """为行星生成卫星系统"""
        # 转换单位
        pmass = planetary_mass * EARTH_MASSES_PER_SOLAR_MASS
        prad = plan_radius / KM_PER_AU
        
        # 计算卫星系统参数
        maxdist = math.sqrt(pmass) / 200  # 最大卫星距离（AU）
        mindist = prad * system.random.range(2.5, 10)
        lastrad = mindist
        
        # 计算最大卫星数量
        maxcount = int(math.sqrt(pmass * 10 + 5)) + 1
        count = 0
        
        # 初始化卫星系统
        head = last = None
        pmass *= system.random.range(0.01, 0.2)
        maxcount = int(system.random.range(maxcount / 10, maxcount))
        maxdist *= system.random.range(0.5, 1.5)
        
        # 生成卫星
        while pmass > 0.001 and count < maxcount and lastrad < maxdist:
            # 计算质量范围
            maxfac = math.sqrt((lastrad - prad) / maxdist) / 8
            massmin = 1e17 / EARTH_MASS_IN_GRAMS
            massmax = system.random.range(pmass / 1e6, pmass * maxfac)
            mmin = math.pow(massmin, 1.0 / 4)
            mmax = math.pow(massmax, 1.0 / 4)
            mass = math.pow(system.random.range(mmin, mmax), 4)
            
            # 计算轨道距离
            dist = math.sqrt(mass) * 50000 / KM_PER_AU
            
            if not (mass > massmin):
                continue
                
            count += 1
            
            # 创建卫星
            moon = Planet()
            moon.mass = mass / EARTH_MASSES_PER_SOLAR_MASS
            moon.a = system.random.range(lastrad, lastrad * 1.3)
            lastrad = moon.a + dist
            moon.e = system.random.eccentricity()
            moon.first_moon = None
            
            pmass -= mass * 2
            
            # 添加到卫星系统
            if last is not None:
                last.next_planet = moon
            else:
                head = moon
            last = moon
            
        return head 
﻿using System;

namespace Accrete
{
    /// <summary>
    ///     The class that simulates the accretation
    /// </summary>
    public class Accretation
    {
        // 设置初始条件：初始化尘埃云和行星系统
        public static void SetInitialConditions(SolarSystem system, Double inner_limit_of_dust, Double outer_limit_of_dust)
        {
            // 创建第一个尘埃带
            system.dust_head = new Dust();
            system.planet_head = null;
            // 设置尘埃带边界和状态
            system.dust_head.next_band = null;
            system.dust_head.outer_edge = outer_limit_of_dust;
            system.dust_head.inner_edge = inner_limit_of_dust;
            system.dust_head.dust_present = true;
            system.dust_head.gas_present = true;
            system.dust_left = true;
            system.cloud_eccentricity = 0.2;
        }

        // 计算恒星的尘埃限制范围
        public static Double StellarDustLimit(SolarSystem system, Double stellar_mass_ratio)
        {
            return 200.0 * Math.Pow(stellar_mass_ratio, 1.0 / 3.0);
        }

        // 计算最内层行星的可能位置
        public static Double InnermostPlanet(Double stellar_mass_ratio)
        {
            return 0.3 * Math.Pow(stellar_mass_ratio, 1.0 / 3.0);
        }

        // 计算最外层行星的可能位置
        public static Double OutermostPlanet(Double stellar_mass_ratio)
        {
            return 50.0 * Math.Pow(stellar_mass_ratio, 1.0 / 3.0);
        }

        // 计算行星的内影响范围（考虑离心率和质量）
        public static Double InnerEffectLimit(SolarSystem system, Double a, Double e, Double mass)
        {
            return a * (1.0 - e) * (1.0 - mass) / (1.0 + system.cloud_eccentricity);
        }

        // 计算行星的外影响范围
        public static Double OuterEffectLimit(SolarSystem system, Double a, Double e, Double mass)
        {
            return a * (1.0 + e) * (1.0 + system.reduced_mass) / (1.0 - system.cloud_eccentricity);
        }

        // 检查指定范围内是否有可用尘埃
        public static Boolean DustAvailable(SolarSystem system, Double inside_range, Double outside_range)
        {
            // 遍历尘埃带，检查是否有尘埃存在
            Dust current_dust_band = system.dust_head;
            while ((current_dust_band != null) && (current_dust_band.outer_edge < inside_range))
                current_dust_band = current_dust_band.next_band;
            Boolean dust_here = current_dust_band != null && current_dust_band.dust_present;
            while ((current_dust_band != null) && (current_dust_band.inner_edge < outside_range))
            {
                dust_here = dust_here || current_dust_band.dust_present;
                current_dust_band = current_dust_band.next_band;
            }
            return dust_here;
        }

        // 更新尘埃带状态
        public static void UpdateDustLanes(ref SolarSystem system, Double min, Double max, Double mass, Double crit_mass, Double body_inner_bound, Double body_outer_bound)
        {
            // 根据行星形成情况更新尘埃带
            system.dust_left = false;
            Boolean gas = !(mass > crit_mass);
            Dust node1 = system.dust_head;
            Dust node2;
            while (node1 != null)
            {
                if ((node1.inner_edge < min) && (node1.outer_edge > max))
                {
                    node2 = new Dust
                    {
                        inner_edge = min,
                        outer_edge = max
                    };
                    if (node1.gas_present)
                        node2.gas_present = gas;
                    node2.dust_present = false;
                    Dust node3 = new Dust
                    {
                        inner_edge = max,
                        outer_edge = node1.outer_edge,
                        gas_present = node1.gas_present,
                        dust_present = node1.dust_present,
                        next_band = node1.next_band
                    };
                    node1.next_band = node2;
                    node2.next_band = node3;
                    node1.outer_edge = min;
                    node1 = node3.next_band;
                }
                else if ((node1.inner_edge < max) && (node1.outer_edge > max))
                {
                    node2 = new Dust
                    {
                        next_band = node1.next_band,
                        dust_present = node1.dust_present,
                        gas_present = node1.gas_present,
                        outer_edge = node1.outer_edge,
                        inner_edge = max
                    };
                    node1.next_band = node2;
                    node1.outer_edge = max;
                    if (node1.gas_present)
                        node1.gas_present = gas;
                    node1.dust_present = false;
                    node1 = node2.next_band;
                }
                else if ((node1.inner_edge < min) && (node1.outer_edge > min))
                {
                    node2 = new Dust
                    {
                        next_band = node1.next_band,
                        dust_present = false
                    };
                    if (node1.gas_present)
                        node2.gas_present = gas;
                    node2.outer_edge = node1.outer_edge;
                    node2.inner_edge = min;
                    node1.next_band = node2;
                    node1.outer_edge = min;
                    node1 = node2.next_band;
                }
                else if ((node1.inner_edge >= min) && (node1.outer_edge <= max))
                {
                    if (node1.gas_present)
                        node1.gas_present = gas;
                    node1.dust_present = false;
                    node1 = node1.next_band;
                }
                else if ((node1.outer_edge < min) || (node1.inner_edge > max))
                    node1 = node1.next_band;
            }
            node1 = system.dust_head;
            while (node1 != null)
            {
                if (node1.dust_present && (node1.outer_edge >= body_inner_bound) && (node1.inner_edge <= body_outer_bound))
                    system.dust_left = true;
                node2 = node1.next_band;
                if (node1.dust_present == node2?.dust_present && (node1.gas_present == node2.gas_present))
                {
                    node1.outer_edge = node2.outer_edge;
                    node1.next_band = node2.next_band;
                    node2 = null;
                }
                node1 = node1.next_band;
            }
        }

        // 收集尘埃物质
        public static Double CollectDust(ref SolarSystem system, Double last_mass, Double a, Double e, Double crit_mass, Dust dust_band)
        {
            Double mass_density;
            Double temp = last_mass / (1.0 + last_mass);
            system.reduced_mass = Math.Pow(temp, 1.0 / 4.0);
            system.r_inner = InnerEffectLimit(system, a, e, system.reduced_mass);
            system.r_outer = OuterEffectLimit(system, a, e, system.reduced_mass);
            if (system.r_inner < 0.0)
                system.r_inner = 0.0;
            if (dust_band == null)
                return 0.0;
            Double temp_density = !dust_band.dust_present ? 0.0 : system.dust_density;
            if ((last_mass < crit_mass) || !dust_band.gas_present)
                mass_density = temp_density;
            else
                mass_density = Constants.K * temp_density / (1.0 + Math.Sqrt(crit_mass / last_mass) * (Constants.K - 1.0));
            if ((dust_band.outer_edge <= system.r_inner) || (dust_band.inner_edge >= system.r_outer))
                return CollectDust(ref system, last_mass, a, e, crit_mass, dust_band.next_band);
                
            Double bandwidth = system.r_outer - system.r_inner;
            Double temp1 = system.r_outer - dust_band.outer_edge;
            if (temp1 < 0.0)
                temp1 = 0.0;
            Double width = bandwidth - temp1;
            Double temp2 = dust_band.inner_edge - system.r_inner;
            if (temp2 < 0.0)
                temp2 = 0.0;
            width = width - temp2;
            temp = 4.0 * Math.PI * Math.Pow(a, 2.0) * system.reduced_mass * (1.0 - e * (temp1 - temp2) / bandwidth);
            Double volume = temp * width;
            return volume * mass_density + CollectDust(ref system, last_mass, a, e, crit_mass, dust_band.next_band);
        }

        /*--------------------------------------------------------------------------*/
        /*   Orbital radius is in AU, eccentricity is unitless, and the stellar     */
        /*  luminosity ratio is with respect to the sun.  The value returned is the */
        /*  mass at which the planet begins to accrete gas as well as dust, and is  */
        /*  in units of solar masses.                                               */
        /*--------------------------------------------------------------------------*/
        // 计算临界质量（开始吸积气体的质量阈值）
        public static Double CriticalLimit(Double orbital_radius, Double eccentricity, Double stellar_luminosity_ratio)
        {
            Double perihelion_dist = orbital_radius - orbital_radius * eccentricity;
            Double temp = perihelion_dist * Math.Sqrt(stellar_luminosity_ratio);
            return Constants.B * Math.Pow(temp, -0.75);
        }

        // 执行尘埃吸积过程
        public static void AccreteDust(ref SolarSystem system, ref Double seed_mass, Double a, Double e, Double crit_mass, Double body_inner_bound, Double body_outer_bound)
        {
            Double temp_mass;
            Double new_mass = seed_mass;
            // 循环吸积尘埃直到质量变化很小
            do
            {
                temp_mass = new_mass;
                new_mass = CollectDust(ref system, new_mass, a, e, crit_mass, system.dust_head);
            } while (!(new_mass - temp_mass < 0.0001 * temp_mass));
            seed_mass = seed_mass + new_mass;
            UpdateDustLanes(ref system, system.r_inner, system.r_outer, seed_mass, crit_mass, body_inner_bound, body_outer_bound);
        }

        // 处理星子合并
        // 处理行星之间的碰撞和合并
        public static void CoalescePlanetesimals(ref SolarSystem system, Double a, Double e, Double mass, Double crit_mass, Double stellar_luminosity_ratio, Double body_inner_bound, Double body_outer_bound)
        {
            Boolean coalesced = false;
            Planet node1 = system.planet_head;
            Planet node2 = null;
            Planet node3 = null;
            while (node1 != null)
            {
                node2 = node1;
                Double temp = node1.a - a;
                Double dist1;
                Double dist2;
                if (temp > 0.0)
                {
                    dist1 = a * (1.0 + e) * (1.0 + system.reduced_mass) - a;
                    /* x aphelion   */
                    system.reduced_mass = Math.Pow(node1.mass / (1.0 + node1.mass), 1.0 / 4.0);
                    dist2 = node1.a
                            - node1.a * (1.0 - node1.e) * (1.0 - system.reduced_mass);
                }
                else
                {
                    dist1 = a - a * (1.0 - e) * (1.0 - system.reduced_mass);
                    /* x perihelion */
                    system.reduced_mass = Math.Pow(node1.mass / (1.0 + node1.mass), 1.0 / 4.0);
                    dist2 = node1.a * (1.0 + node1.e) * (1.0 + system.reduced_mass)
                            - node1.a;
                }
                if ((Math.Abs(temp) <= Math.Abs(dist1)) || (Math.Abs(temp) <= Math.Abs(dist2)))
                {
                    system.Callback("Collision between two planetesimals!\n");
                    Double a3 = (node1.mass + mass) / (node1.mass / node1.a + mass / a);
                    temp = node1.mass * Math.Sqrt(node1.a) * Math.Sqrt(1.0 - Math.Pow(node1.e, 2.0));
                    temp = temp + mass * Math.Sqrt(a) * Math.Sqrt(Math.Sqrt(1.0 - Math.Pow(e, 2.0)));
                    temp = temp / ((node1.mass + mass) * Math.Sqrt(a3));
                    temp = 1.0 - Math.Pow(temp, 2.0);
                    if ((temp < 0.0) || (temp >= 1.0))
                        temp = 0.0;
                    e = Math.Sqrt(temp);
                    temp = node1.mass + mass;
                    AccreteDust(ref system, ref temp, a3, e, stellar_luminosity_ratio, body_inner_bound, body_outer_bound);
                    node1.a = a3;
                    node1.e = e;
                    node1.mass = temp;
                    node1 = null;
                    coalesced = true;
                }
                else
                    node1 = node1.next_planet;
            }
            if (coalesced) return;
            node3 = new Planet
            {
                a = a,
                e = e,
                gas_giant = mass >= crit_mass,
                mass = mass
            };
            system.bodies.Add(node3);
            if (system.planet_head == null)
            {
                system.planet_head = node3;
                node3.next_planet = null;
            }
            else
            {
                node1 = system.planet_head;
                if (a < node1.a)
                {
                    node3.next_planet = node1;
                    system.planet_head = node3;
                }
                else if (system.planet_head.next_planet == null)
                {
                    system.planet_head.next_planet = node3;
                    node3.next_planet = null;
                }
                else
                {
                    while ((node1 != null) && (node1.a < a))
                    {
                        node2 = node1;
                        node1 = node1.next_planet;
                    }
                    node3.next_planet = node1;
                    node2.next_planet = node3;
                }
            }
        }

        // 分配行星质量
        public static Planet DistributePlanetaryMasses(ref SolarSystem system, Double stellar_mass_ratio, Double stellar_luminosity_ratio, Double inner_dust, Double outer_dust)
        // 主要行星形成过程
        // 1. 设置初始条件
        // 2. 循环注入原行星
        // 3. 执行吸积和合并
        {
            SetInitialConditions(system, inner_dust, outer_dust);
            Double planetesimal_inner_bound = InnermostPlanet(stellar_mass_ratio);
            Double planetesimal_outer_bound = OutermostPlanet(stellar_mass_ratio);
            while (system.dust_left)
            {
                Double a = system.random.Range(planetesimal_inner_bound, planetesimal_outer_bound);
                Double e = system.random.Eccentricity();
                Double mass = Constants.PROTOPLANET_MASS;
                if (system.verbose)
                    system.Callback("Checking " + a + " AU.\n");
                if (!DustAvailable(system, InnerEffectLimit(system, a, e, mass), OuterEffectLimit(system, a, e, mass)))
                {
                    if (system.verbose)
                        system.Callback(".. failed.\n");
                    continue;
                }
                system.Callback(".. Injecting protoplanet.\n");
                system.dust_density = Constants.DUST_DENSITY_COEFF * Math.Sqrt(stellar_mass_ratio) * Math.Exp(-Constants.ALPHA * Math.Pow(a, 1.0 / Constants.N));
                Double crit_mass = CriticalLimit(a, e, stellar_luminosity_ratio);
                AccreteDust(ref system, ref mass, a, e, crit_mass, planetesimal_inner_bound, planetesimal_outer_bound);
                if ((mass != 0.0) && (mass != Constants.PROTOPLANET_MASS))
                    CoalescePlanetesimals(ref system, a, e, mass, crit_mass, stellar_luminosity_ratio, planetesimal_inner_bound, planetesimal_outer_bound);
                else system.Callback(".. failed due to large neighbor.\n");
            }
            return system.planet_head;
        }

        // 分配卫星质量
        public static Planet DistributeMoonMasses(ref SolarSystem system, Double planetary_mass, Double plan_radius)
        // 为行星生成卫星系统
        {
            Planet last;
            Double pmass = planetary_mass * Constants.SUN_MASS_IN_EARTH_MASSES;
            Double prad = plan_radius / Constants.KM_PER_AU;
            Double maxdist = Math.Sqrt(pmass) / 200;  /* max. moon distance in AU */
            Double mindist = prad * system.random.Range(2.5, 10);
            Double lastrad = mindist;
            Int32 maxcount = (Int32)Math.Sqrt(pmass * 10 + 5) + 1;
            Int32 count = 0;
            Planet head = last = null;
            pmass *= system.random.Range(0.01, 0.2);
            maxcount = (Int32)system.random.Range(maxcount / 10, maxcount);
            maxdist *= system.random.Range(0.5, 1.5);
            while (pmass > 0.001 && count < maxcount && lastrad < maxdist)
            {
                Double maxfac = Math.Sqrt((lastrad - prad) / maxdist) / 8;
                Double massmin = 1e17 / Constants.EARTH_MASS_IN_GRAMS;
                Double massmax = system.random.Range(pmass / 1e6, pmass * maxfac);
                Double mmin = Math.Pow(massmin, 1.0 / 4);
                Double mmax = Math.Pow(massmax, 1.0 / 4);
                Double mass = Math.Pow(system.random.Range(mmin, mmax), 4);
                Double dist = Math.Sqrt(mass) * 50000 / Constants.KM_PER_AU;
                if (!(mass > massmin)) continue;
                count++;
                Planet moon = new Planet
                {
                    mass = mass / Constants.SUN_MASS_IN_EARTH_MASSES,
                    a = system.random.Range(lastrad, lastrad * 1.3)
                };
                lastrad = moon.a + dist;
                moon.e = system.random.Eccentricity();
                moon.first_moon = null;
                pmass -= mass * 2;
                if (last != null)
                    last.next_planet = moon;
                else
                    head = moon;
                last = moon;
            }
            return head;
        }
    }
}
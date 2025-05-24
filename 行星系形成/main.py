import numpy as np
from celestial_body import DustCloud, CelestialBody, G, AU, EvolutionTime

M_sun = 1

dustcloud = DustCloud(Sigma_0=0.0001913, p=1.5)

def generate_bodies():
    num_bodies = 20
    bodies = []

    r_min, r_max = 0.25, 20.0
    r_samples = np.linspace(r_min, r_max, 1000)
    weights = DustCloud.surface_density(r_samples)
    weights /= weights.sum()

    np.random.seed(42)
    sampled_r = np.random.choice(r_samples, size=num_bodies, p=weights)

    theta_samples = np.random.uniform(0, 2*np.pi, size=num_bodies)
    local_sigma = DustCloud.surface_density(sampled_r)

    area = 0.01  # (0.1 AU)^2，单位面积内吸积量

    mass_samples = np.random.uniform(0.5, 1.5, size=num_bodies) * local_sigma * area  # 质量范围

    density1 = np.random.uniform(6e6, 1.1e7, size=num_bodies/2)
    density2 = np.random.uniform(8.4e6, 3.4e7, size=num_bodies - num_bodies/2)
    density_samples = np.concatenate((density1, density2))

    radius_samples = (3*mass_samples / (4*np.pi*density_samples))**(1/3)  # 假设恒定密度，估算半径

    mean = np.array([0, 0])
    cov = np.array([[0.15, 0], 
                    [0, 0.15]])
    delta_v = np.random.multivariate_normal(mean, cov, size=num_bodies)  # 添加小扰动

    for r, num in zip(sampled_r, range(num_bodies)):
        pos = np.array([r * np.cos(theta_samples[num]), r * np.sin(theta_samples[num])])
        v_kep = np.sqrt(G * M_sun / r)  # 开普勒速度
        v = v_kep * np.array([np.cos(theta_samples[num]), np.sin(theta_samples[num])])+ delta_v[num]  # 添加扰动
        body = CelestialBody(num, mass = mass_samples[num], position = pos, velocity = v, radius = radius_samples[num])
        bodies.append(body)







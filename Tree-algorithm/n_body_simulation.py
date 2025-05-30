import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.spatial import KDTree  # 用于快速邻居搜索

# 参数设置
M_sun = 1.989e30          # 太阳质量 (kg)
r_sun = np.array([0, 0])  # 太阳位置 (原点)
N_init = 100                # 初始粒子数
G = 6.674e-11              # 引力常数
dt = 3e3                   # 时间步长
steps = 1000               # 模拟步数
softening = 1e3            # 软化长度
center = np.array([0.0, 0.0])  # 旋转中心

# 初始化粒子
np.random.seed(42)
mass = np.random.uniform(1e21, 1e22, N_init)
r = np.random.uniform(0.9, 1.1, (N_init,)) * 1.5e11
theta = np.random.uniform(0, 2 * np.pi, (N_init,))
pos = np.zeros((N_init, 2))
pos[:, 0] = r * np.cos(theta)
pos[:, 1] = r * np.sin(theta)

vel = np.zeros_like(pos)
for i in range(N_init):
    r = pos[i] - r_sun
    r_norm = np.linalg.norm(r)
    if r_norm > 0:
        tangent_dir = np.array([-r[1], r[0]]) / r_norm  # 切向方向
        vel[i] = np.sqrt(G * M_sun / r_norm) * tangent_dir

def compute_force(pos, mass, G, softening):
    force = np.zeros_like(pos)
    # 太阳引力
    r_to_sun = pos - r_sun
    dist_to_sun = np.sqrt(np.sum(r_to_sun**2, axis=1))
    force -= G * M_sun * mass[:, np.newaxis] * r_to_sun / dist_to_sun[:, np.newaxis]**3
    
    # 粒子间引力
    for i in range(len(pos)):
        for j in range(len(pos)):
            if i != j:
                r = pos[j] - pos[i]
                dist = np.sqrt(np.sum(r**2) + softening**2)
                force[i] += G * mass[i] * mass[j] * r / dist**3
    return force

def compute_radii(mass):
    """假设密度均匀，半径 R ∝ m^(1/3)"""
    density = 2.5e3  # 假设密度 5000 kg/m³（类似岩石）
    return (3 * mass / (4 * np.pi * density)) ** (1/3)

def handle_collisions(pos, vel, mass):
    if len(mass) < 2:
        return pos, vel, mass
    
    radii = compute_radii(mass)
    tree = KDTree(pos)
    to_remove = set()
    new_pos, new_vel, new_mass = [], [], []
    
    # 预计算所有粒子对的临界距离
    for i in range(len(mass)):
        if i in to_remove:
            continue
        merged = False
        neighbors = tree.query_ball_point(pos[i], r=2 * np.max(radii))  # 初步筛选
        
        for j in neighbors:
            if j <= i or j in to_remove:
                continue
                
            # 计算质量依赖的融合距离
            m_ratio = (mass[i] + mass[j]) / 1e20  # 质量尺度归一化
            d_critical = 2 * (radii[i] + radii[j]) * (1 + m_ratio)  # β=2
            
            # 判定是否融合
            r_ij = np.linalg.norm(pos[i] - pos[j])
            if r_ij < d_critical:
                # 合并粒子（动量守恒）
                total_mass = mass[i] + mass[j]
                new_pos.append((mass[i]*pos[i] + mass[j]*pos[j]) / total_mass)
                new_vel.append((mass[i]*vel[i] + mass[j]*vel[j]) / total_mass)
                new_mass.append(total_mass)
                to_remove.update([i, j])
                merged = True
                break
                
        if not merged:
            new_pos.append(pos[i])
            new_vel.append(vel[i])
            new_mass.append(mass[i])
    
    return np.array(new_pos), np.array(new_vel), np.array(new_mass)

# Verlet积分（适应动态粒子数）
def verlet_integrate(pos, vel, force, mass, dt):
    pos_new = pos + vel * dt + 0.5 * force / mass[:, np.newaxis] * dt**2
    vel_new = vel + 0.5 * force / mass[:, np.newaxis] * dt
    return pos_new, vel_new

# 可视化
fig, ax = plt.subplots(figsize=(8, 8))
scat = ax.scatter(pos[:, 0], pos[:, 1], s=mass/1e21, c=mass, cmap='plasma', alpha=0.7)
ax.set_xlim(-3e11, 3e11)
ax.set_ylim(-3e11, 3e11)
title = ax.set_title(f"N={N_init}, Mass={np.sum(mass):.2e} kg")

def update(frame):
    global pos, vel, mass, force
    force = compute_force(pos, mass, G, softening)
    pos, vel = verlet_integrate(pos, vel, force, mass, dt)
    
    # 处理碰撞并更新粒子
    pos, vel, mass = handle_collisions(pos, vel, mass)
    
    # 更新绘图
    scat.set_offsets(pos)
    scat.set_sizes(mass/1e21)
    scat.set_array(mass)
    title.set_text(f"N={len(mass)}, Mass={np.sum(mass):.2e} kg")
    return scat, title

# 初始力计算
force = compute_force(pos, mass, G, softening)

# 运行动画
ani = FuncAnimation(fig, update, frames=steps, interval=10, blit=False)
plt.colorbar(scat, label='Mass (kg)')
plt.show()

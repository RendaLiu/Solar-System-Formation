import numpy as np
from scipy.spatial import KDTree  # 用于快速构建树
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class TreeNode:
    def __init__(self, bounds, depth=0, max_depth=10, theta=0.5):
        self.bounds = bounds  # [xmin, xmax, ymin, ymax]
        self.children = []
        self.center_of_mass = np.zeros(2)
        self.total_mass = 0
        self.particle = None  # (pos, mass, index)
        self.depth = depth
        self.max_depth = max_depth
        self.theta = theta

    def is_leaf(self):
        return len(self.children) == 0

    def insert(self, pos, mass, index):
        if self.is_leaf():
            if self.particle is None:
                self.particle = (pos, mass, index)
                self.center_of_mass = pos
                self.total_mass = mass
            else:
                self._split()
                self.insert(pos, mass, index)
        else:
            # 更新质心和总质量
            self.center_of_mass = (self.center_of_mass * self.total_mass + pos * mass) / (self.total_mass + mass)
            self.total_mass += mass
            for child in self.children:
                if child._in_bounds(pos):
                    child.insert(pos, mass, index)

    def _split(self):
        xmin, xmax, ymin, ymax = self.bounds
        xmid, ymid = (xmin + xmax) / 2, (ymin + ymax) / 2
        self.children = [
            TreeNode([xmin, xmid, ymin, ymid], self.depth + 1, self.max_depth, self.theta),
            TreeNode([xmid, xmax, ymin, ymid], self.depth + 1, self.max_depth, self.theta),
            TreeNode([xmin, xmid, ymid, ymax], self.depth + 1, self.max_depth, self.theta),
            TreeNode([xmid, xmax, ymid, ymax], self.depth + 1, self.max_depth, self.theta)
        ]
        # 将当前粒子插入子节点
        old_pos, old_mass, old_idx = self.particle
        self.particle = None
        for child in self.children:
            if child._in_bounds(old_pos):
                child.insert(old_pos, old_mass, old_idx)

    def _in_bounds(self, pos):
        x, y = pos
        xmin, xmax, ymin, ymax = self.bounds
        return (xmin <= x <= xmax) and (ymin <= y <= ymax)

    def compute_force_on(self, target_pos, target_mass, target_idx, G):
        force = np.zeros(2)
        if self.is_leaf():
            if self.particle is not None and self.particle[2] != target_idx:  # 排除自身
                r = self.particle[0] - target_pos
                dist = np.sqrt(r[0]**2 + r[1]**2)
                force += G * target_mass * self.particle[1] * r / dist**3
        else:
            r = self.center_of_mass - target_pos
            dist = np.sqrt(r[0]**2 + r[1]**2)
            l = self.bounds[1] - self.bounds[0]  # 节点尺寸
            if l / dist < self.theta:
                force += G * target_mass * self.total_mass * r / dist**3
            else:
                for child in self.children:
                    force += child.compute_force_on(target_pos, target_mass, target_idx, G)
        return force

def build_tree(pos, mass, bounds, theta=0.5):
    """构建Barnes-Hut树"""
    root = TreeNode(bounds, theta=theta)
    for i in range(len(pos)):
        root.insert(pos[i], mass[i], i)
    return root

def compute_force_tree(pos, mass, G, theta=0.5, tidal_k=0.0):
    x_min, x_max = np.min(pos[:, 0]), np.max(pos[:, 0])
    y_min, y_max = np.min(pos[:, 1]), np.max(pos[:, 1])
    bounds = [x_min, x_max, y_min, y_max]
    root = TreeNode(bounds, theta=theta)
    
    # 插入所有粒子并记录索引
    for i in range(len(pos)):
        root.insert(pos[i], mass[i], i)
    
    # 计算每个粒子的受力
    forces = np.zeros_like(pos)
    for i in range(len(pos)):
        forces[i] = root.compute_force_on(pos[i], mass[i], i, G)
        forces[i, 0] += tidal_k * pos[i, 0] * mass[i]
    return forces

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

# Verlet积分（与之前相同）
def verlet_integrate(pos, vel, force, mass, dt, Ly=None):
    pos_new = pos + vel * dt + 0.5 * force / mass[:, np.newaxis] * dt**2
    vel_new = vel + 0.5 * force / mass[:, np.newaxis] * dt
    # y方向周期边界
    if Ly is not None:
        pos_new[:, 1] = (pos_new[:, 1] + Ly / 2) % Ly - Ly / 2
    return pos_new, vel_new

# 初始化粒子
N = 1000
pos = np.random.uniform(-1e9, 1e9, (N, 2))
mass = np.random.uniform(2e21, 1e22, N)
vel = np.random.uniform(-2e2, 2e2, (N, 2))
G = 6.674e-11
dt = 1e3
steps = 1000
Ly = 2e9
M_sun = 1.989e30
AU = 1.496e11
tidal_k = - G * M_sun / AU ** 3

# 可视化
fig, ax = plt.subplots(figsize=(8, 8))
scat = ax.scatter(pos[:, 0], pos[:, 1], s=mass/1e21, c=mass, cmap='plasma', alpha=0.7)
ax.set_xlim(-2e9, 2e9)
ax.set_ylim(-1e9, 1e9)
title = ax.set_title(f"N={N}, Mass={np.sum(mass):.2e} kg")

def update(frame):
    global pos, vel, mass
    force = compute_force_tree(pos, mass, G, theta=0.5, tidal_k=tidal_k)
    pos, vel = verlet_integrate(pos, vel, force, mass, dt, Ly)
    # 处理碰撞并更新粒子
    pos, vel, mass = handle_collisions(pos, vel, mass)
    
    # 更新绘图
    scat.set_offsets(pos)
    scat.set_sizes(mass / 1e21)
    scat.set_array(mass)
    title.set_text(f"N={len(mass)}, Mass={np.sum(mass):.2e} kg")
    return scat, title

ani = FuncAnimation(fig, update, frames=steps, interval=10, blit=False)
plt.colorbar(scat, label='Mass (kg)')
plt.show()

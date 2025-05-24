import numpy as np
AU = 1.495979e+8
EvolutionTime = 10000  # years
SOLAR_MASS = 1
G = 3946.56  # AU^3 / (M_s * yr^2)
''' 
-----------------------------README--------------------------------
一开始AI推荐了DustCloud类自己按照既定函数演化，我推导了这个函数，发现函数的形式比较复杂，于是用以下的logistics函数来近似
def logistics(t):
    return np.exp(-(t/EvolutionTime -1.58)) / (np.exp((t/EvolutionTime-1.58)) + np.exp(-(t/EvolutionTime-1.58))) 
但是后来放弃使用这个思路了，目前的DustCloud与CelestialBody整体是质量守恒的

DustCloud中所有质量损失都表现为整体密度的损失，CelestialBody类的update_mass_velo函数采用动量守恒来更新速度

'''


class DustCloud:
    def __init__(self, Sigma_0=0.0001913, p=1.5, hr=0.05, decay_timescale=np.inf):
        self.Sigma_0 = Sigma_0
        self.p = p
        self.hr = hr
        self.tau_decay = decay_timescale  # e.g. 1e6 yr

    def surface_density(self, r):

        return self.Sigma_0 * (r / 1.0)**(-self.p)
    
    def get_sigma(self, x, y, t=0.0): ## 初始密度函数
        r = np.sqrt(x**2 + y**2)
        sigma_r = self.Sigma_0 * (r / 1.0)**(-self.p) * (r>0.25)
        return sigma_r

    def velocity_field(self, x, y): ## 初始速度函数
        r = np.sqrt(x**2 + y**2) 
        if np.sqrt(x**2 + y**2) > 0: 
            pass
        else:
            return 0.0, 0.0
        h = self.hr * r
        v_theta = np.sqrt(G / r) * (1 - (h / r)**2)**0.25 *(r>0.25) ## M = 1 solar mass
        vx = -v_theta * y / r
        vy = v_theta * x / r
        return vx, vy
    
    def reduce_density(self, x, y, dm):
        r = np.sqrt(x**2 + y**2)
        C = 4*np.pi * (np.sqrt(30) - 0.5) 
        self.Sigma_0 = max(0, self.Sigma_0 - dm / C)
    


def accretion_efficiency(v_rel, v_esc):
    if v_rel < 0.5*v_esc:
        return 0.4 # 高效吸积
    elif v_rel < v_esc:
        return 0.2*(v_esc - v_rel)/(0.5*v_esc)
    else:
        return 0.0 # 无法捕获
    

class CelestialBody(DustCloud):
    def __init__(self, num, mass, position, velocity, radius):
        self.num = num # 原行星编号
        self.mass = mass
        self.pos = np.array(position, dtype=np.float64) 
        self.vel = np.array(velocity, dtype=np.float64) 
        self.radius = radius
    def update_mass_velo(self, dt):
        """
        根据星云密度更新质量
        这里并没有复现引力聚焦效应，当v_rel非常小时，dm/dt非常大，可能会导致数值不稳定，需要寻找进一步资料（目前我没有发现可以直接应用的公式）——lrd
        :param density: 密度 (M_s/AU^2)
        """
        if self.mass > 0:
            Sigma = DustCloud.get_sigma(self, self.pos[0], self.pos[1])
            v_rel = np.linalg.norm(self.vel - DustCloud.velocity_field(self.pos[0], self.pos[1]))  # 相对速度
            v_esc = np.sqrt(2 * G * self.mass / self.radius)  # 逃逸速度
            eta = accretion_efficiency(v_rel, v_esc)
            area = np.pi * self.radius**2  # 假设为球形
            dm = Sigma * area * v_rel * eta * dt

            self.vel = (self.mass * self.vel + dm * DustCloud.velocity_field(self.pos[0], self.pos[1])) / (self.mass + dm)  # 按照动量守恒更新速度
            self.vel += (G / np.linalg.norm(self.pos)**2) * dt  # 更新速度
            self.mass += dm


    def update_position(self, dt):
        """
        根据当前速度更新位置
        :param dt: 时间步长（yr）
        """
        if self.mass > 0:
            self.pos += self.vel * dt
            if np.linalg.norm(self.pos) < 0.01: 
                self.mass = 0 # 如果小于0.01AU，认为已经被太阳吞噬
                self.vel = 0

    def collision(self):
        pass    

    def __str__(self):
        return f"{self.name}: pos={self.position}, vel={self.velocity}"
    


'''
以下为没有完成的grid方法（弃用状态），担心计算量直接爆炸没有使用

class DustCloud:
    def __init__(self, x_grid, y_grid, initial_density_fn, velocity_field_fn=None):
        self.x = x_grid
        self.y = y_grid
        self.X, self.Y = np.meshgrid(x_grid, y_grid)
        
        # 初始密度函数 ρ(x, y, t=0)
        self.density = initial_density_fn(self.X, self.Y)
        
        # 可选：初始速度场 v(x, y)
        if velocity_field_fn is not None:
            self.vx, self.vy = velocity_field_fn(self.X, self.Y)
        else:
            self.vx = np.zeros_like(self.X)
            self.vy = np.zeros_like(self.Y)
        
        self.time = 0.0

    def evolve(self, dt, diffusion_coef=0.0):
        """
        基于简化连续介质方程，更新密度状态。
        支持平移+扩散形式。
        """
        # 速度平移（半拉格朗日法近似）
        shifted_density = np.roll(self.density, shift=-int(self.vx[0,0] * dt), axis=1)
        shifted_density = np.roll(shifted_density, shift=-int(self.vy[0,0] * dt), axis=0)
        
        # 可选扩散项（稳定化）
        if diffusion_coef > 0:
            from scipy.ndimage import gaussian_filter
            shifted_density = gaussian_filter(shifted_density, sigma=diffusion_coef)
        
        self.density = shifted_density
        self.time += dt

    def get_density(self):
        return self.density

    def plot(self):
        import matplotlib.pyplot as plt
        plt.contourf(self.X, self.Y, self.density, levels=50, cmap='plasma')
        plt.colorbar(label='Dust Density')
        plt.title(f"Dust Cloud at t={self.time:.2f}")
        plt.xlabel('x (AU)')
        plt.ylabel('y (AU)')
        plt.axis('equal')
        plt.show()

    def make_density_fn(Sigma_0=0.0001913, p=1.5):
        def density_fn(x, y):
            r = np.sqrt(x**2 + y**2)
            return Sigma_0 * (r / 1.0)**(-p)
        return density_fn

    def make_velocity_field_fn(M_star=1.0):
        def velocity_fn(x, y):
            r = np.sqrt(x**2 + y**2) + 1e-5  # 避免除零
            h = 0.05 * r
            v_phi = np.sqrt(M_star / r) * (1 - (h/r)**2)**0.25
            vx = -v_phi * y / r
            vy =  v_phi * x / r
            return vx, vy
        return velocity_fn

'''
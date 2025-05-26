"""
行星类，用于存储行星的物理和轨道参数
"""

class Planet:
    def __init__(self):
        # 轨道参数
        self.a = 0.0  # 半长轴（AU）
        self.e = 0.0  # 离心率
        self.mass = 0.0  # 质量（太阳质量）
        self.orbit_zone = 0  # 轨道区域
        
        # 物理参数
        self.radius = 0.0  # 半径（km）
        self.density = 0.0  # 密度（g/cm³）
        self.escape_velocity = 0.0  # 逃逸速度（km/s）
        self.molecule_weight = 0.0  # 可保留的最小分子量
        self.rms_velocity = 0.0 
        self.volatile_gas_inventory = 0.0
        self.greenhouse_effect = False
        # 环境参数
        self.surface_accel = 0.0  # 表面加速度（cm/s²）
        self.surface_grav = 0.0  # 表面重力（地球重力）
        self.boil_point = 0.0  # 水的沸点（℃）
        self.surface_pressure = 0.0  # 表面压力（大气压）
        self.surface_temp = 0.0  # 表面温度（℃）
        self.hydrosphere = 0.0  # 水圈比例
        self.cloud_cover = 0.0  # 云层覆盖比例
        self.ice_cover = 0.0  # 冰层覆盖比例
        
        # 其他参数
        self.axial_tilt = 0.0  # 轴向倾斜（度）
        self.albedo = 0.0  # 反照率
        self.orbital_period = 0.0  # 轨道周期（天）
        self.day = 0.0  # 日长（小时）
        
        # 系统参数
        self.gas_giant = False  # 是否为气态巨行星
        self.resonant_period = False  # 是否处于轨道共振

        #卫星链表 后续管理用到时再修改
        self.bodies_orbiting = []
        self.next_planet = None  # 下一个行星
        self.first_moon = None  # 第一个卫星 
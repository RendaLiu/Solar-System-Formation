"""
物理常数和模型参数定义
"""

# 天文单位转换常数
SOLAR_MASS_IN_GRAMS = 1.989e33
EARTH_MASS_IN_GRAMS = 5.977e27
EARTH_MASSES_PER_SOLAR_MASS = 332775.64
KM_PER_AU = 1.496e8
CM_PER_KM = 1e5
CM_PER_METER = 100
CM_PER_AU = 1.495978707e13  # 天文单位对应的厘米数

# 物理常数
GRAV_CONSTANT = 6.672e-8  # 引力常数，单位：dyne cm2/gram2
MOLAR_GAS_CONST = 8314.41  # 气体常数，单位：g*m2/(sec2*K*mol)
EARTH_ACCELERATION = 981.0  # 地球表面重力加速度，单位：cm/sec2
EARTH_EXOSPHERE_TEMP = 1273.0  # 地球外大气层温度，单位：开尔文
EARTH_AXIAL_TILT = 23.45  # 地球自转轴倾角，单位：度
DAYS_IN_A_YEAR = 365.25  # 一年中的天数
SECONDS_PER_HOUR = 3600  # 一小时的秒数
MILLIBARS_PER_BAR = 1000.0  # 1巴等于多少毫巴
KELVIN_CELCIUS_DIFFERENCE = 273.0  # 开尔文温度与摄氏温度的差值

# 行星形成相关参数
PROTOPLANET_MASS = 1.0e-15  # 初始行星质量（太阳质量）
DUST_DENSITY_COEFF = 2.0e-3  # 尘埃密度系数
ALPHA = 5.0  # 尘埃分布参数
N = 3.0  # 尘埃分布指数
K = 50.0  # 气体吸积系数
B = 1.2e-5  # 临界质量系数
ECCENTRICITY_COEFF = 0.077  # 离心率系数

# 地球相关参数
EARTH_RADIUS = 6.378e6  # 地球半径，单位：厘米
EARTH_RADIUS_IN_KM = 6378.0  # 地球半径，单位：千米
EARTH_EFFECTIVE_TEMP = 255.0  # 地球有效温度，单位：开尔文
EARTH_ALBEDO = 0.39  # 地球反照率
EARTH_WATER_MASS_PER_AREA = 3.83e15  # 地球单位面积水质量，单位：克/平方千米
EARTH_SURF_PRES_IN_MILLIBARS = 1000.0  # 地球表面压力，单位：毫巴
EARTH_CONVECTION_FACTOR = 0.43  # 地球对流因子

# 行星表面参数
FREEZING_POINT_OF_WATER = 273.0  # 水的冰点，单位：开尔文
GAS_RETENTION_THRESHOLD = 5.0  # 气体保留阈值
GAS_GIANT_ALBEDO = 0.5  # 气态巨行星反照率
CLOUD_ALBEDO = 0.52  # 云层反照率
AIRLESS_ROCKY_ALBEDO = 0.07  # 无大气岩石行星反照率
ROCKY_ALBEDO = 0.15  # 岩石行星反照率
WATER_ALBEDO = 0.04  # 水面反照率
AIRLESS_ICE_ALBEDO = 0.5  # 无大气冰行星反照率
ICE_ALBEDO = 0.7  # 冰行星反照率
CLOUD_COVERAGE_FACTOR = 1.839e-8  # 云层覆盖因子，单位：Km2/kg

# 环境效应参数
GREENHOUSE_EFFECT_CONST = 0.93  # 温室效应常数
J = 1.46e-19  # 日长计算常数，单位：cm2/sec2 g

# 行星形成计算参数
A1_20 = 6.485e12  # Kothari半径计算参数
A2_20 = 4.0032e-8  # Kothari半径计算参数
BETA_20 = 5.71e12  # Kothari半径计算参数
Q1_36 = 1.258e19  # 云层覆盖计算参数，单位：克
Q2_36 = 0.0698  # 云层覆盖计算参数，单位：1/开尔文

# 分子量阈值
MOLECULAR_HYDROGEN = 2.0  # H2
HELIUM = 4.0  # He
METHANE = 16.0  # CH4
AMMONIA = 17.0  # NH3
WATER_VAPOR = 18.0  # H2O
NEON = 20.2  # Ne
MOLECULAR_NITROGEN = 28.0  # N2
CARBON_MONOXIDE = 28.0  # CO
NITRIC_OXIDE = 30.0  # NO
MOLECULAR_OXYGEN = 32.0  # O2
HYDROGEN_SULPHIDE = 34.1  # H2S
ARGON = 39.9  # Ar
CARBON_DIOXIDE = 44.0  # CO2
NITROUS_OXIDE = 44.0  # N2O
NITROGEN_DIOXIDE = 46.0  # NO2
OZONE = 48.0  # O3
SULPHUR_DIOXIDE = 64.1  # SO2
SULPHUR_TRIOXIDE = 80.1  # SO3
KRYPTON = 83.8  # Kr
XENON = 131.3  # Xe

# 其他常数
INCREDIBLY_LARGE_NUMBER = 9.9999e37  # 用于特殊计算的极大数

# 新增的原子量常数
ATOMIC_HYDROGEN = 1.0    # H
ATOMIC_NITROGEN = 14.0   # N
ATOMIC_OXYGEN = 16.0     # O

#
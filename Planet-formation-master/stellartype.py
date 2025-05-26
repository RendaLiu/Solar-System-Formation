"""
恒星类型类，用于定义和获取恒星类型
"""

# 先定义类
class StellarType:
    def __init__(self, star_class="", temp=0.0, balmer="", lines="", 
                 mass=0.0, size=0.0, density=0.0, lum=0.0, star_age=0.0):
        self.star_class = star_class  # 恒星光谱类型
        self.temp = temp              # 表面温度
        self.balmer = balmer          # 巴尔默线强度
        self.lines = lines            # 光谱线特征
        self.mass = mass              # 质量（太阳质量）
        self.size = size              # 大小（太阳半径）
        self.density = density        # 密度（太阳密度）
        self.lum = lum                # 光度（太阳光度）
        self.star_age = star_age      # 寿命（年）

    @staticmethod
    def get_stellar_type_mass(mass):
        """根据质量获取恒星类型"""
        for star_type in reversed(StellarType.BUILTIN):
            if mass <= star_type.mass:
                return star_type
        return StellarType.BUILTIN[-1]  # 默认返回最后一个类型

    @staticmethod
    def get_stellar_type_temp(temperature):
        """根据温度获取恒星类型"""
        for star_type in reversed(StellarType.BUILTIN):
            if temperature <= star_type.temp:
                return star_type
        return StellarType.BUILTIN[-1]  # 默认返回最后一个类型

# 在类定义后创建BUILTIN列表
StellarType.BUILTIN = [
    StellarType("O0", 1e10, "weak", "He+ O-II He-II", 40, 17.8, 0.01, 405000, 1e6),
    StellarType("B0", 30000, "medium", "He", 18, 7.4, 0.1, 13000, 11e6),
    StellarType("A0", 12000, "strong", "", 3.5, 2.5, 0.3, 80, 440e6),
    StellarType("F0", 7500, "medium", "", 1.7, 1.4, 1.0, 6.4, 3e9),
    StellarType("G0", 6000, "weak", "Ca++ Fe++", 1.1, 1.0, 1.4, 1.4, 8e9),
    StellarType("K0", 5000, "v. weak", "Ca++ Fe++", 0.8, 0.8, 1.8, 0.46, 17e9),
    StellarType("M0", 3500, "v. weak", "Ca++ TiO2", 0.5, 0.6, 2.5, 0.08, 56e9),
    StellarType("D0", 1500, "none", "", 0, 0, 2.5, 0.00, 56e9)
]

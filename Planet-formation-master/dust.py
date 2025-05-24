"""
尘埃带类，用于模拟行星形成过程中的尘埃云结构
"""

class Dust:
    def __init__(self):
        self.inner_edge = 0.0  # 内边界（AU）
        self.outer_edge = 0.0  # 外边界（AU）
        self.dust_present = True  # 是否存在尘埃
        self.gas_present = True   # 是否存在气体
        self.next_band = None     # 下一个尘埃带 
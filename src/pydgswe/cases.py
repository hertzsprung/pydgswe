from . import BoundaryConditions

def three_humps(xx):
    zz = 0.0

    if xx >= 22.0 and xx < 25.0:
        zz = (0.05)*xx - 1.1
    elif xx >=25.0 and xx <= 28.0:
        zz = (-0.05)*xx + 1.4
    elif xx > 8.0 and xx < 12.0:
        zz = 0.2 - (0.05*(xx-10)**2)
    elif xx > 39.0 and xx < 46.5:
        zz = 0.3
    else:
        zz = 0.0
    
    return 10.0*zz

class LakeAtRest:
    def __init__(self):
        self.xmin = 0.0
        self.xmax = 50.0
        self.bcs = BoundaryConditions()
        self.manning = 0.0
        self.simulation_time = 0.5
        self.Bed_Data = three_humps

    def Init_Conds_h(self, zz, xx):
        return 2.0 - zz

    def Init_Conds_q(self, zz, xx):
        return 0.0

class BuildingOvertopping:
    def __init__(self):
        self.xmin = 0.0
        self.xmax = 50.0
        self.bcs = BoundaryConditions()
        self.manning = 0.02
        self.simulation_time = 10.0
        self.Bed_Data = three_humps

    def Init_Conds_h(self, zz, xx):
        return 6 - zz if xx <= 25 else 0

    def Init_Conds_q(self, zz, xx):
        return 0.0

class SheetFlow:
    def __init__(self):
        self.xmin = 0.0
        self.xmax = 1000.0
        self.bcs = BoundaryConditions()
        self.bcs.reflect_up = -1
        self.bcs.reflect_down = -1
        self.manning = 0.035
        self.simulation_time = 600.0
        self.zs = [291.263908, 291.18313652271297, 290.8903204624281, 290.47561599699134, 289.995009942357, 289.4431302663146, 288.82279961126136, 288.1374964177651, 287.3031160487284, 286.1799920732785, 284.758743390925, 283.32592827613877, 281.8881232097677, 280.44311468575836, 279.0193793524901, 277.9049990848545, 277.62594722961046, 273.5178066125639, 276.632811909643, 281.242813456406, 282.20686260695163, 283.10937521226873, 284.29406801761803, 285.9109340880744, 287.8456278186377, 289.62905870985406, 290.13156192017976, 289.9793708587651, 289.8000032881645, 289.61999462791044, 289.4400016055666, 289.2934420746389, 289.4078067731277, 288.65937841886114, 287.55220061503303, 286.3421937579219, 285.82907237347456, 287.0581357115665, 287.76969897501425, 288.2331243078665, 288.44967626009225, 288.6637574047714, 289.3725128699639, 290.19218498512737, 291.16000444268775, 292.5106190691995, 294.3075102550974, 296.70344444749765, 299.5212554509419, 302.654373, 302.654373]

    def Bed_Data(self, xx):
        i = int(xx / 20.0)
        return self.zs[i]
        
    def Init_Conds_h(self, zz, xx):
        return 1.5e-3
#        if xx > 800.0:
#            return max(0.0, 290.0-zz)
#        else:
#            return 0.0

    def Init_Conds_q(self, zz, xx):
        return 0.0

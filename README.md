# 玩具矩阵库v1


这个分支是用来测试"自动复数"的
自动复数就是统计矩阵里面每个元素的虚部，然后判断这个矩阵需不需要复数运算
在STM32G4(Arm Cortex-M4)平台上，一次浮点除法运算需要消耗14个时钟周期，因此有必要判断是不是在执行高斯消元的时候开启复数除法。
# 玩具矩阵库v1


这个分支是用来测试"自动复数"的

自动复数就是统计矩阵里面每个元素的虚部，然后判断这个矩阵需不需要复数运算

在STM32G4(Arm Cortex-M4)平台上，一次浮点除法运算需要消耗14个时钟周期，因此有必要判断是不是在执行高斯消元的时候开启复数除法。

This branch is aimed to test "auto complex switch", this function is dedicated to do statistics on the imaginary part of each entry that composing a matrix, then give a judgement on whether the matrix needs complex arithmetic. 

On STM32G4(a SoC based on Arm Cortex-M4 Architecture), a single float-point number divide operation costs 14 clock cycle, much enough to give a judge on enabling complex arithmetic during Gauss Elimination.

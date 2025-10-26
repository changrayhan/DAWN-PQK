#ifndef FAST_INVERSION_H
#define FAST_INVERSION_H

#include "Poly.h"
#include <optional>

/**
 * DAWN 论文 Algorithm 1: FastInversion 函数
 * 计算 f^{-1} mod (x^{main_n/2} + 1, 2) 在 Z_2 上的逆元
 * 
 * 该函数实现了 DAWN 论文中的快速逆元算法，利用恒等式 
 * x^{2^l} + 1 = (x+1)^{2^l} 在 Z_2 中成立的性质。
 * 
 * @param f 输入多项式（整数多项式）
 * @param main_n 主环度数（必须是2的幂），算法计算 f^{-1} mod (x^{main_n/2} + 1, 2)
 * @return std::optional<Poly> 如果可逆则返回逆元，否则返回 std::nullopt
 * @throws std::invalid_argument 如果 main_n 不是2的幂
 * 
 * 算法步骤：
 * 1. 检查 f mod (x+1, 2) != 0，否则不可逆
 * 2. 计算 k = (f + 1) / (x + 1)（在 Z_2 中的多项式长除法）
 * 3. 初始化 f2 = 1
 * 4. 循环 i = 0 到 log2(main_n/2) - 1：
 *    - b = (k * f2) mod (x^{2^i} + 1, 2)
 *    - k = (k + f * b) mod (x^{main_n/2} + 1, 2)
 *    - k = k / (x^{2^i} + 1)（利用多项式除法）
 *    - f2 = f2 + (x^{2^i} + 1) * b
 * 5. 返回 f2
 */
std::optional<Poly> fast_inversion(const Poly& f, int main_n);

#endif // FAST_INVERSION_H

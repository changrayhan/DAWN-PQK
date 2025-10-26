#include "FastInversion.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

// 前向声明辅助函数
Poly divide_by_x_plus_1(const Poly& p, int mod);
Poly divide_by_xm_plus_1(const Poly& k, int m, int mod);

/**
 * DAWN 论文 Algorithm 1: FastInversion 函数
 * 计算 f^{-1} mod (x^{main_n} + 1, 2) 在 Z_2 上的逆元
 * 
 * 注意：注释中写 "main_n/2" 是误导性的，实际计算的是 mod (x^{main_n} + 1, 2)
 * 当传入 main_n = n/4 时，实际上计算 mod (x^{n/4} + 1, 2)，这是正确的
 * 
 * @param f 输入多项式（整数多项式）
 * @param main_n 目标子环度数（必须是2的幂），算法计算 f^{-1} mod (x^{main_n} + 1, 2)
 * @return std::optional<Poly> 如果可逆则返回逆元，否则返回 std::nullopt
 */
std::optional<Poly> fast_inversion(const Poly& f, int main_n) {
    // 参数验证
    if (main_n <= 0 || (main_n & (main_n - 1)) != 0) {
        throw std::invalid_argument("main_n 必须是2的幂");
    }
    
    int target_n = main_n; // 目标子环度数
    int log_target = static_cast<int>(std::log2(target_n)); // 循环次数
    
    // 步骤1: 检查 f mod (x+1, 2) != 0，否则不可逆
    Poly f_mod_x_plus_1 = f.mod_xm_plus_1(1, 2);
    if (f_mod_x_plus_1.is_zero()) {
        return std::nullopt; // 不可逆
    }
    
    // 步骤2: 计算 k = (f + 1) / (x + 1)（在 Z_2 中的多项式长除法）
    Poly f_plus_1 = f + Poly::one_poly(f.get_n(), f.get_q());
    Poly k = divide_by_x_plus_1(f_plus_1, 2);
    
    // 步骤3: 初始化 f2 = 1 (in Z2[x]/(x^{main_n}+1))
    Poly f2 = Poly::one_poly(f.get_n(), 2);
    
    // 步骤4: 循环 i = 0 到 log2(target_n) - 1
    for (int i = 0; i < log_target; ++i) {
        int m = 1 << i; // m = 2^i
        
        // b = (k * f2) mod (x^{2^i} + 1, 2)
        Poly b = (k * f2).mod_xm_plus_1_keep_degree(m, 2);
        
        // k = (k + f * b) mod (x^{target_n} + 1, 2) [保持度数 main_n]
        Poly f_mod_z2 = f.mod_xm_plus_1_keep_degree(f.get_n(), 2);
        // 将 f_mod_z2 转换为 Z_2 多项式
        Poly f_z2(f.get_n(), 2);
        for (int j = 0; j < f.get_n(); ++j) {
            f_z2.set_coeff(j, f_mod_z2.get_coeff(j) % 2);
        }
        Poly fb = (f_z2 * b).mod_xm_plus_1_keep_degree(target_n, 2);
        k = (k + fb).mod_xm_plus_1_keep_degree(target_n, 2);
        
        // k = k / (x^{2^i} + 1)（利用多项式除法）
        k = divide_by_xm_plus_1(k, m, 2);
        
        // f2 = f2 + (x^{2^i} + 1) * b
        Poly xm_plus_1 = Poly::xm_plus_1(f.get_n(), 2, m);
        f2 = (f2 + xm_plus_1 * b).mod_xm_plus_1_keep_degree(target_n, 2);
    }
    
    // 步骤5: 返回 f2 mod (x^{target_n} + 1, 2)
    return f2.mod_xm_plus_1(target_n, 2);
}

/**
 * 辅助函数：多项式除以 (x + 1) 在 Z_2 上
 * 在 Z_2[x] 中，p(x) / (x+1) 的系数可通过递推计算：
 * q[0] = p[0]
 * q[i] = p[i] + q[i-1] for i >= 1
 */
Poly divide_by_x_plus_1(const Poly& p, int mod) {
    int n = p.get_n();
    Poly q(n, mod); // 使用传入的模数
    
    // 将输入多项式的系数转换为指定模数
    q.set_coeff(0, p.get_coeff(0) % mod);
    for (int i = 1; i < n; ++i) {
        int val = (p.get_coeff(i) + q.get_coeff(i-1)) % mod;
        if (val < 0) val += mod;
        q.set_coeff(i, val);
    }
    
    return q;
}

/**
 * 辅助函数：多项式除以 (x^m + 1) 在 Z_2 上
 * 前提：k 是 (x^m + 1) 的倍数
 * 在 Z2 中，x^m ≡ 1，所以 k(x) = (x^m + 1) * q(x) ⇒
 * k[i] = q[i] + q[i - m] （下标模 n）
 * 因此 q[i] = k[i] + q[i - m]
 */
Poly divide_by_xm_plus_1(const Poly& k, int m, int mod) {
    // 前提：k 是 (x^m + 1) 的倍数
    int n = k.get_n();
    Poly q(n, 2); // 强制使用 q=2，所有运算在 Z_2 中
    
    for (int i = 0; i < n; ++i) {
        int val = k.get_coeff(i);
        if (i >= m) {
            val ^= q.get_coeff(i - m); // Z2 中加法 = XOR
        }
        q.set_coeff(i, val);
    }
    
    return q;
}
#include "Poly.h"
#include "FastInversion.h"
#include "TernarySampler.h"
#include <vector>
#include <cstdint>
#include <algorithm>
#include <climits>
#include <iostream>

/**
 * DAWN 论文 Algorithm 2: SimpleDecoding 函数
 * 
 * 功能：从编码的多项式中恢复原始消息
 * 
 * @param c_prime 在 Z_q 中的多项式（长度 n）
 * @param m_prime 在 Z_2 中的多项式（长度 n/2）
 * @param f2 FastInversion 的结果（在 Z_2 中，长度 n/4）
 * @param params DAWN 参数
 * @return 恢复的原始消息（n/4 位，以字节向量表示）
 */
std::vector<uint8_t> simple_decoding(const Poly& c_prime, const Poly& m_prime, const Poly& f2, const DAWNParams& params) {
    int n = params.n;
    int q = params.q;
    int n_over_4 = n / 4;
    int n_over_2 = n / 2;
    
    // 调试：打印输入参数
    std::cout << "SimpleDecoding 调试信息:" << std::endl;
    std::cout << "n = " << n << ", q = " << q << ", n/4 = " << n_over_4 << std::endl;
    
    // 步骤 1: 计算 e_prime = c_prime mod (x^{n/4} + 1, 2)
    // 即取前 n/4 系数异或后 n/4 系数
    Poly e_prime = c_prime.mod_xm_plus_1(n_over_4, 2);
    
    // 调试：打印 e_prime 的前几个系数
    std::cout << "e_prime 前8个系数: ";
    for (int i = 0; i < std::min(8, n_over_4); ++i) {
        std::cout << (int)e_prime.get_coeff(i) << " ";
    }
    std::cout << std::endl;
    
    // 步骤 2: 检查 e_prime 是否全零
    bool e_prime_is_zero = true;
    for (int i = 0; i < n_over_4; ++i) {
        if (e_prime.get_coeff(i) != 0) {
            e_prime_is_zero = false;
            break;
        }
    }
    
    std::cout << "e_prime 是否全零: " << (e_prime_is_zero ? "是" : "否") << std::endl;
    
    Poly m(n_over_2, 2); // 在 Z_2 中
    
    if (e_prime_is_zero) {
        // 步骤 2: 若 e_prime 全零，则 m = m_prime mod x^{n/4}
        for (int i = 0; i < n_over_4; ++i) {
            m.set_coeff(i, m_prime.get_coeff(i));
        }
    } else {
        // 步骤 3: 否则进行错误修正
        int min_v = INT_MAX;
        int best_j = 0;
        
        // 对 j = 0,1,2,3，计算 idx = j * n/4
        for (int j = 0; j < 4; ++j) {
            int idx = j * n_over_4;
            int16_t coeff = c_prime.get_coeff(idx);
            
            // 计算 v = min(|c_prime[idx] + q|, |c_prime[idx] - q|)
            int v_plus_q = std::abs(coeff + q);
            int v_minus_q = std::abs(coeff - q);
            int v = std::min(v_plus_q, v_minus_q);
            
            if (v < min_v) {
                min_v = v;
                best_j = j;
            }
        }
        
        // 找到最小 v 对应的 j，设 i = j * n/4
        int i = best_j * n_over_4;
        
        // 修正：m = (m_prime + x^i * f2 * e_prime) mod (x^{n/2} + 1, 2)
        // 使用 Poly 类的乘法运算计算 f2 * e_prime（在 Z_2 中，模 x^{n/4} + 1）
        Poly f2_e_prime = (f2 * e_prime).mod_xm_plus_1(n_over_4, 2);
        
        // 计算 x^i * f2 * e_prime（在 Z_2 中，模 x^{n/2} + 1）
        Poly x_i_f2_e_prime(n_over_2, 2);
        for (int k = 0; k < n_over_4; ++k) {
            int target_idx = (i + k) % n_over_2;
            x_i_f2_e_prime.set_coeff(target_idx, f2_e_prime.get_coeff(k));
        }
        
        // 计算 m = m_prime + x^i * f2 * e_prime（在 Z_2 中）
        for (int k = 0; k < n_over_2; ++k) {
            int sum = m_prime.get_coeff(k) ^ x_i_f2_e_prime.get_coeff(k);
            m.set_coeff(k, sum & 1);
        }
    }
    
    // 步骤 4: 将 n/4 位打包为 std::vector<uint8_t>（每字节 8 位，低位在前）
    std::vector<uint8_t> result;
    int bits_per_byte = 8;
    int bytes_needed = (n_over_4 + bits_per_byte - 1) / bits_per_byte;
    result.resize(bytes_needed, 0);
    
    for (int i = 0; i < n_over_4; ++i) {
        int byte_idx = i / bits_per_byte;
        int bit_idx = i % bits_per_byte;
        if (m.get_coeff(i) & 1) {
            result[byte_idx] |= (1 << bit_idx);
        }
    }
    
    return result;
}

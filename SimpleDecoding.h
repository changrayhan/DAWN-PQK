#ifndef SIMPLE_DECODING_H
#define SIMPLE_DECODING_H

#include "Poly.h"
#include "TernarySampler.h"
#include <vector>
#include <cstdint>

/**
 * DAWN 论文 Algorithm 2: SimpleDecoding 函数
 * 
 * 功能：从编码的多项式中恢复原始消息
 * 
 * @param c_prime 在 Z_q 中的多项式（长度 n），**注意：系数范围应为 [-q, q]，未中心化**
 * @param m_prime 在 Z_2 中的多项式（长度 n/2）
 * @param f2 FastInversion 的结果（在 Z_2 中，长度 n/4）
 * @param params DAWN 参数
 * @return 恢复的原始消息（n/4 位，以字节向量表示）
 * 
 * 算法步骤：
 * 1. 计算 e_prime = c_prime mod (x^{n/4} + 1, 2)
 * 2. 若 e_prime 全零，则 m = m_prime mod x^{n/4}
 * 3. 否则进行错误修正：
 *    - 对 j = 0,1,2,3，计算 idx = j * n/4
 *    - 计算 v = min(|c_prime[idx] + q|, |c_prime[idx] - q|)
 *    - 找到最小 v 对应的 j，设 i = j * n/4
 *    - 修正：m = (m_prime + x^i * f2 * e_prime) mod (x^{n/2} + 1, 2)
 *    - 取前 n/4 位作为消息
 * 4. 将 n/4 位打包为 std::vector<uint8_t>（每字节 8 位，低位在前）
 * 
 * 重要说明：
 * - c_prime 的系数范围必须为 [-q, q]，不能是中心化的 [-q/2, q/2]
 * - 这是因为错误修正需要检测 c_prime[idx] 是否接近 ±q
 */
std::vector<uint8_t> simple_decoding(const Poly& c_prime, const Poly& m_prime, const Poly& f2, const DAWNParams& params);

#endif // SIMPLE_DECODING_H

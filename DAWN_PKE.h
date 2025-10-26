#ifndef DAWN_PKE_H
#define DAWN_PKE_H

#include "Poly.h"
#include "FastInversion.h"
#include "TernarySampler.h"
#include "SimpleDecoding.h"
#include <vector>
#include <cstdint>
#include <utility>
#include <optional>

/**
 * DAWN.PKE 公钥加密体制实现
 * 
 * 该类实现了 DAWN 论文中的三个核心算法：
 * - KeyGen (Algorithm 3): 密钥生成
 * - Encrypt (Algorithm 4): 加密
 * - Decrypt (Algorithm 5): 解密
 * 
 * 基于多项式环 R = Z_q[x]/(x^n + 1) 和子环运算
 */
class DAWN_PKE {
private:
    DAWNParams params;  // DAWN 参数
    
    /**
     * 预计算 w = x^{n/4} + 1 多项式
     * 用于加密时的消息编码
     */
    Poly w_poly;
    
    /**
     * 预计算 t = x^{n/2} + 1 多项式
     * 用于解密时的子环运算
     */
    Poly t_poly;
    
    /**
     * 检查多项式 f 在环 R_{x^n+1,q} 和 R_{x^{n/2}+1,2} 中是否均可逆
     * 通过实际计算逆元来验证可逆性
     * @param f 待检查的多项式
     * @return 如果均可逆返回 true，否则返回 false
     */
    bool is_invertible(const Poly& f) const;
    
    /**
     * 将消息字节向量转换为多项式
     * @param message 消息字节向量（长度为 n/4 字节）
     * @return 对应的多项式（系数为 0/1，长度为 n/4，高位补零至 n）
     */
    Poly message_to_poly(const std::vector<uint8_t>& message) const;
    
    /**
     * 从种子字节向量生成两个子种子
     * @param seed 原始种子字节向量
     * @return std::pair<std::vector<uint8_t>, std::vector<uint8_t>> 两个子种子
     */
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> split_seed(const std::vector<uint8_t>& seed) const;

public:
    /**
     * 构造函数
     * @param params DAWN 参数
     */
    explicit DAWN_PKE(const DAWNParams& params);
    
    /**
     * DAWN 论文 Algorithm 3: KeyGen
     * 
     * 密钥生成算法：
     * 1. 循环采样 f, g 直到 f 在 R_{x^n+1,q} 和 R_{x^{n/2}+1,2} 中均可逆
     * 2. 计算 h = g * f^{-1} mod (x^n+1, q)
     * 3. 调用 fast_inversion(f, n/4) 得 f2
     * 
     * @return std::pair<Poly, std::pair<Poly, Poly>> 
     *         (pk, (f, f2)) 其中 pk 是公钥，f 是私钥，f2 是快速逆元
     */
    std::pair<Poly, std::pair<Poly, Poly>> keygen();
    
    /**
     * DAWN 论文 Algorithm 4: Encrypt
     * 
     * 加密算法：
     * 1. 将 message（n/4 位）扩展为多项式 m
     * 2. 用 seed 拆分为 seed_s, seed_e
     * 3. 采样 s ~ T_{n,ks}, e ~ T_{n,ke}
     * 4. 计算 c = h*s + e + w*m（其中 w = x^{n/4}+1）
     * 5. 调用 c.compress(dc)
     * 
     * @param pk 公钥多项式
     * @param message 消息字节向量（长度为 n/4 字节）
     * @param seed 随机种子字节向量
     * @return 密文多项式
     */
    Poly encrypt(const Poly& pk, const std::vector<uint8_t>& message, const std::vector<uint8_t>& seed) const;
    
    /**
     * DAWN 论文 Algorithm 5: Decrypt
     * 
     * 解密算法：
     * 1. 解压缩密文 c
     * 2. 计算 c_prime = t * c * f mod (x^n+1, q)
     * 3. 计算 m_prime = (c_prime * f2) mod (x^{n/2}+1, 2)
     * 4. 将 c_prime 约简到 n/2 度数
     * 5. 调用 simple_decoding(c_prime_for_decoding, m_prime, f2, params)
     * 
     * @param sk 私钥对 (f, f2)
     * @param c 密文多项式
     * @return 解密后的消息字节向量
     */
    std::vector<uint8_t> decrypt(const std::pair<Poly, Poly>& sk, const Poly& c) const;
    
    /**
     * 获取 DAWN 参数
     * @return DAWNParams 参数结构体
     */
    const DAWNParams& get_params() const { return params; }
    
    /**
     * 获取消息长度（以字节为单位）
     * @return 消息长度
     */
    int get_message_length() const { return (params.n / 4 + 7) / 8; }  // 将位转换为字节
    
    // 调试方法
    Poly message_to_poly_debug(const std::vector<uint8_t>& message) const { return message_to_poly(message); }
    
    // 更多调试方法
    Poly get_t_poly() const { return t_poly; }
    Poly get_w_poly() const { return w_poly; }
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> split_seed_debug(const std::vector<uint8_t>& seed) const { return split_seed(seed); }
};

#endif // DAWN_PKE_H

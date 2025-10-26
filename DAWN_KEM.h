#ifndef DAWN_KEM_H
#define DAWN_KEM_H

#include "DAWN_PKE.h"
#include "Poly.h"
#include <vector>
#include <cstdint>
#include <utility>
#include <openssl/evp.h>
#include <openssl/crypto.h>

/**
 * DAWN.KEM 密钥封装机制实现
 * 
 * 该类实现了 DAWN 论文中的 Algorithms 6-8：
 * - KeyGen (Algorithm 6): 密钥生成
 * - Encaps (Algorithm 7): 密钥封装
 * - Decaps (Algorithm 8): 密钥解封装
 * 
 * 基于 DAWN.PKE 公钥加密体制和 SHAKE256 哈希函数
 */
class DAWN_KEM {
public:
    /**
     * 使用 SHAKE256 对公钥进行哈希
     * @param pk 公钥多项式
     * @param out_len 输出长度（字节）
     * @return 哈希结果
     */
    static std::vector<uint8_t> hash_pk(const Poly& pk, size_t out_len);
    
    /**
     * 使用 SHAKE256 对消息进行哈希，生成密钥和随机种子
     * @param m 消息字节向量
     * @param hpk 公钥哈希值
     * @param out_len 输出长度（字节），默认为96（32字节K + 64字节rho）
     * @return std::pair<K, rho> 其中 K 是密钥，rho 是随机种子
     */
    static std::pair<std::vector<uint8_t>, std::vector<uint8_t>> hash_m(
        const std::vector<uint8_t>& m, 
        const std::vector<uint8_t>& hpk, 
        size_t out_len = 96);
    
    /**
     * 恒定时间比较两个多项式是否相等
     * 用于实现隐式拒绝，抵抗 CCA 攻击
     * @param a 第一个多项式
     * @param b 第二个多项式
     * @return 如果相等返回 true，否则返回 false
     */
    static bool constant_time_compare(const Poly& a, const Poly& b);

private:
    DAWN_PKE pke;  // DAWN PKE 实例
    
    /**
     * 使用 SHAKE256 对密钥和密文进行哈希
     * @param K 密钥字节向量
     * @param c 密文多项式
     * @param out_len 输出长度（字节），默认为32
     * @return 哈希结果
     */
    static std::vector<uint8_t> hash_K(
        const std::vector<uint8_t>& K, 
        const Poly& c, 
        size_t out_len = 32);

public:
    /**
     * 私钥结构体
     * 包含 PKE 私钥、公钥哈希值和随机密钥
     */
    struct PrivateKey {
        std::pair<Poly, Poly> pke_sk;  // PKE 私钥 (f, f2)
        Poly pk;                        // 公钥
        std::vector<uint8_t> H_pk;      // 公钥哈希值
        std::vector<uint8_t> K_prime;   // 随机密钥 K'
        
        PrivateKey(std::pair<Poly, Poly> sk, Poly public_key, 
                  std::vector<uint8_t> hpk, std::vector<uint8_t> k_prime)
            : pke_sk(std::move(sk)), pk(std::move(public_key)), 
              H_pk(std::move(hpk)), K_prime(std::move(k_prime)) {}
    };
    
    /**
     * 构造函数
     * @param params DAWN 参数
     */
    explicit DAWN_KEM(const DAWNParams& params);
    
    /**
     * DAWN 论文 Algorithm 6: KeyGen
     * 
     * 密钥生成算法：
     * 1. 调用 DAWN_PKE::keygen() 生成 (pk, (f, f2))
     * 2. 生成随机 K_prime（n/4 位）
     * 3. 计算 H_pk = hash_pk(pk)
     * 
     * @return std::pair<Poly, PrivateKey> (pk, sk)
     */
    std::pair<Poly, PrivateKey> keygen();
    
    /**
     * DAWN 论文 Algorithm 7: Encaps
     * 
     * 密钥封装算法：
     * 1. 随机生成 m（n/4 位）
     * 2. (K, rho) = hash_m(m, H_pk)
     * 3. c = pke.encrypt(pk, m, bytes_to_seed(rho))
     * 4. K = hash_K(K, c)
     * 
     * @param pk 公钥多项式
     * @return std::pair<Poly, std::vector<uint8_t>> (c, K) 密文和密钥
     */
    std::pair<Poly, std::vector<uint8_t>> encaps(const Poly& pk);
    
    /**
     * DAWN 论文 Algorithm 8: Decaps
     * 
     * 密钥解封装算法：
     * 1. m = pke.decrypt(sk.pke_sk, c)
     * 2. (K, rho) = hash_m(m, sk.H_pk)
     * 3. 重新加密 c_recomputed = pke.encrypt(sk.pk, m, bytes_to_seed(rho))
     * 4. 若 c == c_recomputed，返回 hash_K(K, c)；否则返回 hash_K(sk.K_prime, c)
     * 
     * @param sk 私钥
     * @param c 密文多项式
     * @return 解封装的密钥
     */
    std::vector<uint8_t> decaps(const PrivateKey& sk, const Poly& c);
    
    /**
     * 获取 DAWN 参数
     * @return DAWNParams 参数结构体
     */
    const DAWNParams& get_params() const { return pke.get_params(); }
    
    /**
     * 获取消息长度（以字节为单位）
     * @return 消息长度
     */
    int get_message_length() const { return pke.get_message_length(); }
    
    // 调试方法
    const DAWN_PKE& get_pke() const { return pke; }
    static std::vector<uint8_t> hash_K_debug(
        const std::vector<uint8_t>& K, 
        const Poly& c, 
        size_t out_len = 32) { return hash_K(K, c, out_len); }
};

#endif // DAWN_KEM_H

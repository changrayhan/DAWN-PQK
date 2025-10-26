#ifndef TERNARY_SAMPLER_H
#define TERNARY_SAMPLER_H

#include "Poly.h"
#include <vector>
#include <random>
#include <optional>
#include <stdexcept>
#include <string>
#include <cstdint>

/**
 * DAWN 参数管理结构体
 * 包含 DAWN 加密体制的所有关键参数
 */
struct DAWNParams {
    int n;           // 多项式度数（必须是2的幂）
    int q;           // 模数
    int kg;          // 密钥生成参数
    int kf;          // 前向安全参数
    int ks;          // 签名参数
    int ke;          // 加密参数
    int dc;          // 压缩位数
    std::string name; // 参数集名称
    
    /**
     * 构造函数
     * @param n 多项式度数
     * @param q 模数
     * @param kg 密钥生成参数
     * @param kf 前向安全参数
     * @param ks 签名参数
     * @param ke 加密参数
     * @param dc 压缩位数
     * @param name 参数集名称
     */
    DAWNParams(int n, int q, int kg, int kf, int ks, int ke, int dc, const std::string& name)
        : n(n), q(q), kg(kg), kf(kf), ks(ks), ke(ke), dc(dc), name(name) {}
};

/**
 * 三元分布采样器类
 * 用于生成符合三元分布的多项式（系数为 -1, 0, +1）
 */
class TernarySampler {
public:
    /**
     * 三元分布采样函数
     * 生成长度为 n 的多项式，其中 +1 和 -1 各 k 个，其余为 0
     * 
     * @param n 多项式长度
     * @param q 模数
     * @param k +1 和 -1 的个数（各 k 个）
     * @param seed 可选的随机种子，用于确定性采样
     * @return Poly 对象（系数为 -1, 0, +1）
     * @throws std::invalid_argument 如果 2*k > n
     */
    static Poly sample_ternary(int n, int q, int k, std::optional<uint64_t> seed = std::nullopt);
    
    /**
     * 三元分布采样函数（使用字节向量种子）
     * 生成长度为 n 的多项式，其中 +1 和 -1 各 k 个，其余为 0
     * 
     * @param n 多项式长度
     * @param q 模数
     * @param k +1 和 -1 的个数（各 k 个）
     * @param seed 随机种子字节向量
     * @return Poly 对象（系数为 -1, 0, +1）
     * @throws std::invalid_argument 如果 2*k > n
     */
    static Poly sample_ternary(int n, int q, int k, const std::vector<uint8_t>& seed);
    
    /**
     * 使用 DAWN 参数进行三元分布采样
     * 
     * @param params DAWN 参数结构体
     * @param k +1 和 -1 的个数（各 k 个）
     * @param seed 可选的随机种子
     * @return Poly 对象
     */
    static Poly sample_ternary(const DAWNParams& params, int k, std::optional<uint64_t> seed = std::nullopt);
    
    /**
     * 采样 f 多项式（使用 kf 参数）
     * 
     * @param params DAWN 参数结构体
     * @param seed 可选的随机种子
     * @return Poly 对象
     */
    static Poly sample_f(const DAWNParams& params, std::optional<uint64_t> seed = std::nullopt);
    
    /**
     * 采样 f 多项式（使用 kf 参数，字节向量种子）
     * 
     * @param params DAWN 参数结构体
     * @param seed 随机种子字节向量
     * @return Poly 对象
     */
    static Poly sample_f(const DAWNParams& params, const std::vector<uint8_t>& seed);
    
    /**
     * 采样 g 多项式（使用 kg 参数）
     * 
     * @param params DAWN 参数结构体
     * @param seed 可选的随机种子
     * @return Poly 对象
     */
    static Poly sample_g(const DAWNParams& params, std::optional<uint64_t> seed = std::nullopt);
    
    /**
     * 采样 g 多项式（使用 kg 参数，字节向量种子）
     * 
     * @param params DAWN 参数结构体
     * @param seed 随机种子字节向量
     * @return Poly 对象
     */
    static Poly sample_g(const DAWNParams& params, const std::vector<uint8_t>& seed);
    
    /**
     * 采样 s 多项式（使用 ks 参数）
     * 
     * @param params DAWN 参数结构体
     * @param seed 可选的随机种子
     * @return Poly 对象
     */
    static Poly sample_s(const DAWNParams& params, std::optional<uint64_t> seed = std::nullopt);
    
    /**
     * 采样 s 多项式（使用 ks 参数，字节向量种子）
     * 
     * @param params DAWN 参数结构体
     * @param seed 随机种子字节向量
     * @return Poly 对象
     */
    static Poly sample_s(const DAWNParams& params, const std::vector<uint8_t>& seed);
    
    /**
     * 采样 e 多项式（使用 ke 参数）
     * 
     * @param params DAWN 参数结构体
     * @param seed 可选的随机种子
     * @return Poly 对象
     */
    static Poly sample_e(const DAWNParams& params, std::optional<uint64_t> seed = std::nullopt);
    
    /**
     * 采样 e 多项式（使用 ke 参数，字节向量种子）
     * 
     * @param params DAWN 参数结构体
     * @param seed 随机种子字节向量
     * @return Poly 对象
     */
    static Poly sample_e(const DAWNParams& params, const std::vector<uint8_t>& seed);

private:
    /**
     * 使用 Fisher-Yates shuffle 算法选择不重复的索引
     * 
     * @param n 总长度
     * @param k 选择的个数
     * @param rng 随机数生成器
     * @return 选中的索引向量
     */
    static std::vector<int> fisher_yates_sample(int n, int k, std::mt19937_64& rng);
    
    /**
     * 使用 Fisher-Yates shuffle 算法选择不重复的索引（使用字节向量）
     * 
     * @param n 总长度
     * @param k 选择的个数
     * @param random_bytes 随机字节向量
     * @return 选中的索引向量
     */
    static std::vector<int> fisher_yates_sample_with_bytes(int n, int k, const std::vector<uint8_t>& random_bytes);
};

// 预定义的 DAWN 参数集（基于论文 Table 6）
extern const DAWNParams DAWN_ALPHA_512;
extern const DAWNParams DAWN_ALPHA_1024;
extern const DAWNParams DAWN_BETA_512;
extern const DAWNParams DAWN_BETA_1024;

#endif // TERNARY_SAMPLER_H

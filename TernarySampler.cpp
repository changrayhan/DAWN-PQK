#include "TernarySampler.h"
#include <algorithm>
#include <random>
#include <chrono>
#include <openssl/evp.h>
#include <numeric>

// 预定义的 DAWN 参数集（基于论文 Table 6）
// DAWN-α-512: n=512, q=769, kg=160, kf=64, ks=96, ke=160, dc=7
const DAWNParams DAWN_ALPHA_512(512, 769, 160, 64, 96, 160, 7, "DAWN-α-512");
// DAWN-α-1024: n=1024, q=769, kg=256, kf=96, ks=192, ke=256, dc=4
const DAWNParams DAWN_ALPHA_1024(1024, 769, 256, 96, 192, 256, 4, "DAWN-α-1024");
// DAWN-β-512: n=512, q=257, kg=64, kf=32, ks=48, ke=64, dc=2
const DAWNParams DAWN_BETA_512(512, 257, 64, 32, 48, 64, 2, "DAWN-β-512");
// DAWN-β-1024: n=1024, q=257, kg=96, kf=64, ks=96, ke=96, dc=1
const DAWNParams DAWN_BETA_1024(1024, 257, 96, 64, 96, 96, 1, "DAWN-β-1024");

Poly TernarySampler::sample_ternary(int n, int q, int k, std::optional<uint64_t> seed) {
    // 边界条件检查
    if (2 * k > n) {
        throw std::invalid_argument("错误：2*k (" + std::to_string(2 * k) + 
                                  ") 不能大于 n (" + std::to_string(n) + ")");
    }
    
    if (n <= 0 || k < 0 || q <= 0) {
        throw std::invalid_argument("错误：n, k, q 必须为正数");
    }
    
    // 初始化随机数生成器
    std::mt19937_64 rng;
    if (seed.has_value()) {
        rng.seed(seed.value());
    } else {
        // 使用随机设备生成种子
        std::random_device rd;
        rng.seed(rd());
    }
    
    // 创建长度为 n 的零向量
    std::vector<int16_t> coeffs(n, 0);
    
    // 使用 Fisher-Yates shuffle 选择不重复的索引
    std::vector<int> indices = fisher_yates_sample(n, 2 * k, rng);
    
    // 前 k 个设为 +1，后 k 个设为 -1
    for (int i = 0; i < k; ++i) {
        coeffs[indices[i]] = 1;
    }
    for (int i = k; i < 2 * k; ++i) {
        coeffs[indices[i]] = -1;
    }
    
    // 构造并返回 Poly 对象，使用传入的模数 q
    return Poly(coeffs, q, n);
}

Poly TernarySampler::sample_ternary(const DAWNParams& params, int k, std::optional<uint64_t> seed) {
    return sample_ternary(params.n, params.q, k, seed);
}

// 专用采样函数
Poly TernarySampler::sample_f(const DAWNParams& params, std::optional<uint64_t> seed) {
    return sample_ternary(params.n, params.q, params.kf, seed);
}

Poly TernarySampler::sample_g(const DAWNParams& params, std::optional<uint64_t> seed) {
    return sample_ternary(params.n, params.q, params.kg, seed);
}

Poly TernarySampler::sample_s(const DAWNParams& params, std::optional<uint64_t> seed) {
    return sample_ternary(params.n, params.q, params.ks, seed);
}

Poly TernarySampler::sample_ternary(int n, int q, int k, const std::vector<uint8_t>& seed) {
    // 边界条件检查
    if (2 * k > n) {
        throw std::invalid_argument("错误：2*k (" + std::to_string(2 * k) + 
                                  ") 不能大于 n (" + std::to_string(n) + ")");
    }
    
    if (n <= 0 || k < 0 || q <= 0) {
        throw std::invalid_argument("错误：n, k, q 必须为正数");
    }
    
    // 使用 SHAKE256 从种子生成随机数
    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    if (!ctx) {
        throw std::runtime_error("无法创建 EVP_MD_CTX");
    }
    
    const EVP_MD* md = EVP_shake256();
    if (!EVP_DigestInit_ex(ctx, md, nullptr)) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法初始化 SHAKE256");
    }
    
    if (!EVP_DigestUpdate(ctx, seed.data(), seed.size())) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法更新 SHAKE256");
    }
    
    // 生成足够的随机字节用于采样
    // 每个索引需要 log2(n) 位，总共需要 2*k*log2(n) 位
    size_t bytes_needed = (2 * k * 32 + 7) / 8;  // 保守估计，使用32位表示索引
    std::vector<uint8_t> random_bytes(bytes_needed);
    if (!EVP_DigestFinalXOF(ctx, random_bytes.data(), bytes_needed)) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法完成 SHAKE256");
    }
    
    EVP_MD_CTX_free(ctx);
    
    // 创建长度为 n 的零向量
    std::vector<int16_t> coeffs(n, 0);
    
    // 使用 Fisher-Yates shuffle 选择不重复的索引
    std::vector<int> indices = fisher_yates_sample_with_bytes(n, 2 * k, random_bytes);
    
    // 前 k 个设为 +1，后 k 个设为 -1
    for (int i = 0; i < k; ++i) {
        coeffs[indices[i]] = 1;
    }
    for (int i = k; i < 2 * k; ++i) {
        coeffs[indices[i]] = -1;
    }
    
    // 构造并返回 Poly 对象，使用传入的模数 q
    return Poly(coeffs, q, n);
}

Poly TernarySampler::sample_s(const DAWNParams& params, const std::vector<uint8_t>& seed) {
    return sample_ternary(params.n, params.q, params.ks, seed);
}

Poly TernarySampler::sample_e(const DAWNParams& params, const std::vector<uint8_t>& seed) {
    return sample_ternary(params.n, params.q, params.ke, seed);
}

Poly TernarySampler::sample_f(const DAWNParams& params, const std::vector<uint8_t>& seed) {
    return sample_ternary(params.n, params.q, params.kf, seed);
}

Poly TernarySampler::sample_g(const DAWNParams& params, const std::vector<uint8_t>& seed) {
    return sample_ternary(params.n, params.q, params.kg, seed);
}

std::vector<int> TernarySampler::fisher_yates_sample(int n, int k, std::mt19937_64& rng) {
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    
    // Fisher-Yates shuffle 算法
    for (int i = 0; i < k; ++i) {
        std::uniform_int_distribution<int> dist(i, n - 1);
        int j = dist(rng);
        std::swap(indices[i], indices[j]);
    }
    
    // 返回前 k 个元素
    indices.resize(k);
    return indices;
}

std::vector<int> TernarySampler::fisher_yates_sample_with_bytes(int n, int k, const std::vector<uint8_t>& random_bytes) {
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    
    size_t byte_index = 0;
    
    // Fisher-Yates shuffle 算法
    for (int i = 0; i < k; ++i) {
        // 生成 [i, n-1] 范围内的随机数
        int range = n - i;
        int bits_needed = 0;
        int temp = range - 1;
        while (temp > 0) {
            bits_needed++;
            temp >>= 1;
        }
        
        int j;
        do {
            // 从随机字节中提取足够的位
            uint32_t random_value = 0;
            int bits_extracted = 0;
            
            while (bits_extracted < bits_needed && byte_index < random_bytes.size()) {
                random_value |= (static_cast<uint32_t>(random_bytes[byte_index]) << bits_extracted);
                bits_extracted += 8;
                byte_index++;
            }
            
            // 如果字节不够，循环使用现有字节
            if (byte_index >= random_bytes.size()) {
                byte_index = 0;
            }
            
            j = i + (random_value % range);
        } while (j < i || j >= n);
        
        std::swap(indices[i], indices[j]);
    }
    
    // 返回前 k 个元素
    indices.resize(k);
    return indices;
}

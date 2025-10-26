#include "DAWN_PKE.h"
#include <random>
#include <stdexcept>
#include <iostream>
#include <openssl/evp.h>
#include <numeric>

// 构造函数
DAWN_PKE::DAWN_PKE(const DAWNParams& params) : params(params), w_poly(params.n, params.q), t_poly(params.n, params.q) {
    // 预计算 w = x^{n/4} + 1 多项式
    w_poly.set_coeff(0, 1);  // x^0 的系数
    w_poly.set_coeff(params.n / 4, 1);  // x^{n/4} 的系数
    
    // 预计算 t = x^{n/2} + 1 多项式
    t_poly.set_coeff(0, 1);  // x^0 的系数
    t_poly.set_coeff(params.n / 2, 1);  // x^{n/2} 的系数
}

// 检查多项式 f 在环 R_{x^n+1,q} 和 R_{x^{n/2}+1,2} 中是否均可逆
bool DAWN_PKE::is_invertible(const Poly& f) const {
    // 检查在 R_{x^{n/2}+1,2} 中是否可逆（使用 fast_inversion）
    auto f2_opt = fast_inversion(f, params.n / 2);
    if (!f2_opt.has_value()) {
        return false;
    }
    
    // 检查在 R_{x^n+1,q} 中是否可逆（使用 Poly::inv）
    auto f_inv_opt = f.inv();
    if (!f_inv_opt.has_value()) {
        return false;
    }
    
    return true;
}

// 将消息字节向量转换为多项式
Poly DAWN_PKE::message_to_poly(const std::vector<uint8_t>& message) const {
    int n_over_4 = params.n / 4;  // n/4 位
    int message_bytes = (n_over_4 + 7) / 8;  // 转换为字节数
    
    // 检查消息长度
    if (message.size() != static_cast<size_t>(message_bytes)) {
        throw std::invalid_argument("消息长度必须为 " + std::to_string(message_bytes) + " 字节");
    }
    
    // 创建长度为 n 的多项式，前 n/4 个系数来自消息，其余为0
    std::vector<int16_t> coeffs(params.n, 0);
    
    // 将消息的每一位转换为多项式系数
    for (int i = 0; i < n_over_4; ++i) {
        int byte_idx = i / 8;
        int bit_idx = i % 8;
        if (byte_idx < static_cast<int>(message.size())) {
            // 提取第 i 位
            coeffs[i] = (message[byte_idx] >> bit_idx) & 1;
        }
    }
    
    return Poly(coeffs, params.q, params.n);
}

// 从种子字节向量生成两个子种子
std::pair<std::vector<uint8_t>, std::vector<uint8_t>> DAWN_PKE::split_seed(const std::vector<uint8_t>& seed) const {
    // 使用 SHAKE256 扩展种子
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
    
    // 生成64字节输出，分为两个32字节的子种子
    std::vector<uint8_t> output(64);
    if (!EVP_DigestFinalXOF(ctx, output.data(), 64)) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法完成 SHAKE256");
    }
    
    EVP_MD_CTX_free(ctx);
    
    std::vector<uint8_t> seed_s(output.begin(), output.begin() + 32);
    std::vector<uint8_t> seed_e(output.begin() + 32, output.end());
    
    return {seed_s, seed_e};
}

// DAWN 论文 Algorithm 3: KeyGen
std::pair<Poly, std::pair<Poly, Poly>> DAWN_PKE::keygen() {
    Poly f(params.n, params.q);
    Poly g(params.n, params.q);
    
    // 循环采样 f, g 直到 f 在 R_{x^n+1,q} 和 R_{x^{n/2}+1,2} 中均可逆
    int max_attempts = 1000;  // 防止无限循环
    int attempts = 0;
    
    do {
        f = TernarySampler::sample_f(params);
        g = TernarySampler::sample_g(params);
        attempts++;
        
        // 如果常数项为0，强制设为1
        if (f.get_coeff(0) == 0) {
            f.set_coeff(0, 1);
        }
        
        if (attempts > max_attempts) {
            throw std::runtime_error("密钥生成失败：无法找到可逆的 f 多项式");
        }
    } while (!is_invertible(f));
    
    // 计算 h = g * f^{-1} mod (x^n+1, q)
    auto f_inv_opt = f.inv();
    if (!f_inv_opt.has_value()) {
        throw std::runtime_error("密钥生成失败：无法计算 f 的逆元");
    }
    Poly f_inv = f_inv_opt.value();
    Poly h = g * f_inv;
    
    // 调用 fast_inversion(f, n/4) 得 f2
    auto f2_opt = fast_inversion(f, params.n / 4);
    if (!f2_opt.has_value()) {
        throw std::runtime_error("密钥生成失败：无法计算 f 的快速逆元");
    }
    Poly f2 = f2_opt.value();
    
    // 返回 (pk, (f, f2))
    return {h, {f, f2}};
}

// DAWN 论文 Algorithm 4: Encrypt
Poly DAWN_PKE::encrypt(const Poly& pk, const std::vector<uint8_t>& message, const std::vector<uint8_t>& seed) const {
    // 将 message（n/4 位）扩展为多项式 m
    Poly m = message_to_poly(message);
    
    // 用 seed 拆分为 seed_s, seed_e
    auto [seed_s, seed_e] = split_seed(seed);
    
    // 采样 s ~ T_{n,ks}, e ~ T_{n,ke}
    Poly s = TernarySampler::sample_s(params, seed_s);
    Poly e = TernarySampler::sample_e(params, seed_e);
    
    // 调试：打印 e 的前几个系数
    std::cout << "错误多项式 e 前8个系数: ";
    for (int i = 0; i < std::min(8, params.n); ++i) {
        std::cout << (int)e.get_coeff(i) << " ";
    }
    std::cout << std::endl;
    
    // 计算 c = h*s + e + w*m（其中 w = x^{n/4}+1）
    Poly w_m = w_poly * m;
    
    // 调试：打印 w_m 的前几个系数
    std::cout << "w_m 前8个系数: ";
    for (int i = 0; i < std::min(8, params.n); ++i) {
        std::cout << (int)w_m.get_coeff(i) << " ";
    }
    std::cout << std::endl;
    
    Poly c = pk * s + e + w_m;
    
    // 调用 c.compress(dc)
    c.compress(params.dc);
    
    return c;
}

// DAWN 论文 Algorithm 5: Decrypt
std::vector<uint8_t> DAWN_PKE::decrypt(const std::pair<Poly, Poly>& sk, const Poly& c) const {
    const Poly& f = sk.first;
    const Poly& f2 = sk.second;
    
    // 解压缩密文 c
    Poly c_decompressed = c;
    c_decompressed.decompress(params.dc, params.q);
    
    // 计算 c_prime = t * c * f mod (x^n+1, q)
    Poly c_prime_temp = t_poly * c_decompressed * f;
    
    // 关键修复：手动恢复未中心化的系数（范围 [0, q)）
    // 避免 Poly 构造函数自动调用 reduce_centered()
    std::vector<int16_t> c_prime_coeffs = c_prime_temp.get_coeffs();
    for (auto& coeff : c_prime_coeffs) {
        // 将系数从 [-q/2, q/2] 范围恢复到 [0, q) 范围
        if (coeff < 0) {
            coeff += params.q;
        }
        // 确保系数在 [0, q) 范围内
        coeff = coeff % params.q;
    }
    
    // 创建未中心化的 c_prime
    Poly c_prime(c_prime_coeffs, params.q, params.n, true);
    
    // 调试：打印 c_prime 的关键系数
    std::cout << "调试 c_prime 系数:" << std::endl;
    std::cout << "c_prime[0] = " << c_prime.get_coeff(0) << std::endl;
    std::cout << "c_prime[" << params.n/4 << "] = " << c_prime.get_coeff(params.n/4) << std::endl;
    std::cout << "c_prime[" << params.n/2 << "] = " << c_prime.get_coeff(params.n/2) << std::endl;
    std::cout << "c_prime[" << 3*params.n/4 << "] = " << c_prime.get_coeff(3*params.n/4) << std::endl;
    
    // 计算 m_prime = (c_prime * f2) mod (x^{n/2}+1, 2)
    // 先将c_prime的系数模2转换到Z2
    Poly c_prime_z2(params.n, 2);
    for (int i = 0; i < params.n; ++i) {
        int16_t coeff = c_prime.get_coeff(i);
        c_prime_z2.set_coeff(i, (coeff % 2 + 2) % 2); // 确保结果在{0,1}中
    }
    
    // 将 f2 扩展到 n 度
    Poly f2_extended(params.n, 2);
    for (int i = 0; i < params.n / 4; ++i) {
        f2_extended.set_coeff(i, f2.get_coeff(i));
    }
    
    // 使用 Poly 类的内置方法计算 c_prime_z2 * f2_extended mod (x^{n/2}+1, 2)
    Poly m_prime = (c_prime_z2 * f2_extended).mod_xm_plus_1(params.n / 2, 2);
    
    // 调用 simple_decoding，传入未中心化的c_prime
    return simple_decoding(c_prime, m_prime, f2, params);
}


#include "DAWN_KEM.h"
#include <random>
#include <stdexcept>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <numeric>

// 构造函数
DAWN_KEM::DAWN_KEM(const DAWNParams& params) : pke(params) {
    // 初始化 OpenSSL
    OpenSSL_add_all_algorithms();
}

// 使用 SHAKE256 对公钥进行哈希
std::vector<uint8_t> DAWN_KEM::hash_pk(const Poly& pk, size_t out_len) {
    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    if (!ctx) {
        throw std::runtime_error("无法创建 EVP_MD_CTX");
    }
    
    const EVP_MD* md = EVP_shake256();
    if (!EVP_DigestInit_ex(ctx, md, nullptr)) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法初始化 SHAKE256");
    }
    
    // 将多项式系数序列化为字节
    const auto& coeffs = pk.get_coeffs();
    std::vector<uint8_t> pk_bytes;
    pk_bytes.reserve(coeffs.size() * sizeof(int16_t));
    
    for (int16_t coeff : coeffs) {
        // 将 int16_t 转换为字节（小端序）
        pk_bytes.push_back(static_cast<uint8_t>(coeff & 0xFF));
        pk_bytes.push_back(static_cast<uint8_t>((coeff >> 8) & 0xFF));
    }
    
    if (!EVP_DigestUpdate(ctx, pk_bytes.data(), pk_bytes.size())) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法更新 SHAKE256");
    }
    
    std::vector<uint8_t> result(out_len);
    if (!EVP_DigestFinalXOF(ctx, result.data(), out_len)) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法完成 SHAKE256");
    }
    
    EVP_MD_CTX_free(ctx);
    return result;
}

// 使用 SHAKE256 对消息进行哈希，生成密钥和随机种子
std::pair<std::vector<uint8_t>, std::vector<uint8_t>> DAWN_KEM::hash_m(
    const std::vector<uint8_t>& m, 
    const std::vector<uint8_t>& hpk, 
    size_t out_len) {
    
    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    if (!ctx) {
        throw std::runtime_error("无法创建 EVP_MD_CTX");
    }
    
    const EVP_MD* md = EVP_shake256();
    if (!EVP_DigestInit_ex(ctx, md, nullptr)) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法初始化 SHAKE256");
    }
    
    // 输入消息 m
    if (!EVP_DigestUpdate(ctx, m.data(), m.size())) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法更新 SHAKE256 (m)");
    }
    
    // 输入公钥哈希 hpk
    if (!EVP_DigestUpdate(ctx, hpk.data(), hpk.size())) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法更新 SHAKE256 (hpk)");
    }
    
    std::vector<uint8_t> result(out_len);
    if (!EVP_DigestFinalXOF(ctx, result.data(), out_len)) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法完成 SHAKE256");
    }
    
    EVP_MD_CTX_free(ctx);
    
    // 将结果分为两部分：K（32字节）和 rho（剩余字节）
    size_t k_len = 32;  // K 固定为32字节
    size_t rho_len = out_len - k_len;  // rho 为剩余字节
    
    std::vector<uint8_t> K(result.begin(), result.begin() + k_len);
    std::vector<uint8_t> rho(result.begin() + k_len, result.end());
    
    return {K, rho};
}

// 使用 SHAKE256 对密钥和密文进行哈希
std::vector<uint8_t> DAWN_KEM::hash_K(
    const std::vector<uint8_t>& K, 
    const Poly& c, 
    size_t out_len) {
    
    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    if (!ctx) {
        throw std::runtime_error("无法创建 EVP_MD_CTX");
    }
    
    const EVP_MD* md = EVP_shake256();
    if (!EVP_DigestInit_ex(ctx, md, nullptr)) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法初始化 SHAKE256");
    }
    
    // 输入密钥 K
    if (!EVP_DigestUpdate(ctx, K.data(), K.size())) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法更新 SHAKE256 (K)");
    }
    
    // 将密文多项式序列化为字节
    const auto& coeffs = c.get_coeffs();
    std::vector<uint8_t> c_bytes;
    c_bytes.reserve(coeffs.size() * sizeof(int16_t));
    
    for (int16_t coeff : coeffs) {
        // 将 int16_t 转换为字节（小端序）
        c_bytes.push_back(static_cast<uint8_t>(coeff & 0xFF));
        c_bytes.push_back(static_cast<uint8_t>((coeff >> 8) & 0xFF));
    }
    
    if (!EVP_DigestUpdate(ctx, c_bytes.data(), c_bytes.size())) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法更新 SHAKE256 (c)");
    }
    
    std::vector<uint8_t> result(out_len);
    if (!EVP_DigestFinalXOF(ctx, result.data(), out_len)) {
        EVP_MD_CTX_free(ctx);
        throw std::runtime_error("无法完成 SHAKE256");
    }
    
    EVP_MD_CTX_free(ctx);
    return result;
}

// 恒定时间比较两个多项式是否相等
bool DAWN_KEM::constant_time_compare(const Poly& a, const Poly& b) {
    // 检查度数和模数是否相同
    if (a.get_n() != b.get_n() || a.get_q() != b.get_q()) {
        return false;
    }
    
    const auto& coeffs_a = a.get_coeffs();
    const auto& coeffs_b = b.get_coeffs();
    
    // 使用 CRYPTO_memcmp 进行恒定时间比较
    return CRYPTO_memcmp(coeffs_a.data(), coeffs_b.data(), 
                        coeffs_a.size() * sizeof(int16_t)) == 0;
}

// DAWN 论文 Algorithm 6: KeyGen
std::pair<Poly, DAWN_KEM::PrivateKey> DAWN_KEM::keygen() {
    // 1. 调用 DAWN_PKE::keygen() 生成 (pk, (f, f2))
    auto [pk, pke_sk] = pke.keygen();
    
    // 2. 生成随机 K_prime（n/4 位）
    const DAWNParams& params = pke.get_params();
    int message_bytes = pke.get_message_length();
    std::vector<uint8_t> K_prime(message_bytes);
    
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_int_distribution<uint8_t> dist(0, 255);
    
    for (int i = 0; i < message_bytes; ++i) {
        K_prime[i] = dist(rng);
    }
    
    // 3. 计算 H_pk = hash_pk(pk)
    std::vector<uint8_t> H_pk = hash_pk(pk, 32);  // 使用32字节输出
    
    // 构造私钥
    PrivateKey sk(pke_sk, pk, H_pk, K_prime);
    
    return {pk, sk};
}

// DAWN 论文 Algorithm 7: Encaps
std::pair<Poly, std::vector<uint8_t>> DAWN_KEM::encaps(const Poly& pk) {
    const DAWNParams& params = pke.get_params();
    int message_bytes = pke.get_message_length();
    
    // 1. 随机生成 m（n/4 位）
    std::vector<uint8_t> m(message_bytes);
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_int_distribution<uint8_t> dist(0, 255);
    
    for (int i = 0; i < message_bytes; ++i) {
        m[i] = dist(rng);
    }
    
    // 2. (K, rho) = hash_m(m, H_pk)
    // 注意：这里需要 H_pk，但在封装时我们只有 pk
    // 根据论文，我们需要先计算 H_pk = hash_pk(pk)
    std::vector<uint8_t> H_pk = hash_pk(pk, 32);
    auto [K, rho] = hash_m(m, H_pk, 96);  // 输出96字节，分为32字节K和64字节rho
    
    // 3. c = pke.encrypt(pk, m, rho)
    Poly c = pke.encrypt(pk, m, rho);
    
    // 4. K = hash_K(K, c)
    std::vector<uint8_t> final_K = hash_K(K, c, 32);
    
    return {c, final_K};
}

// DAWN 论文 Algorithm 8: Decaps
std::vector<uint8_t> DAWN_KEM::decaps(const PrivateKey& sk, const Poly& c) {
    // 1. m = pke.decrypt(sk.pke_sk, c)
    std::vector<uint8_t> m = pke.decrypt(sk.pke_sk, c);
    
    // 2. (K, rho) = hash_m(m, sk.H_pk)
    auto [K, rho] = hash_m(m, sk.H_pk, 96);  // 输出96字节，分为32字节K和64字节rho
    
    // 3. 重新加密 c_recomputed = pke.encrypt(sk.pk, m, rho)
    // 注意：这里使用 sk.pk 而不是传入的 pk，因为解封装时我们只有私钥
    Poly c_recomputed = pke.encrypt(sk.pk, m, rho);
    
    // 4. 若 c == c_recomputed，返回 hash_K(K, c)；否则返回 hash_K(sk.K_prime, c)
    if (constant_time_compare(c, c_recomputed)) {
        return hash_K(K, c, 32);
    } else {
        return hash_K(sk.K_prime, c, 32);
    }
}

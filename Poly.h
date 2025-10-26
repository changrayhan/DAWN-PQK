#ifndef POLY_H
#define POLY_H

#include <vector>
#include <stdexcept>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <climits>
#include <optional>

/**
 * 多项式类，用于支持 DAWN 加密体制中的环运算
 * 在环 R = Z_q[x]/(x^n + 1) 上进行运算
 */
class Poly {
private:
    std::vector<int16_t> coeffs;  // 多项式系数
    int q;                        // 模数
    int n;                        // 多项式度数（必须是2的幂）

    // 检查 n 是否为 2 的幂
    static bool is_power_of_two(int n);
    
    // 将系数中心化到 [-q/2, q/2]
    void reduce_centered();

public:
    /**
     * 构造函数：从系数向量初始化
     * @param coeffs 系数向量，长度必须等于 n
     * @param q 模数（必须为正）
     * @param n 多项式度数（必须是2的幂）
     * @throws std::invalid_argument 如果参数无效
     */
    Poly(const std::vector<int16_t>& coeffs, int q, int n);
    
    /**
     * 默认构造函数（创建零多项式）
     * @param n 多项式度数
     * @param q 模数
     */
    Poly(int n, int q);
    
    /**
     * 构造函数：从原始系数向量初始化（不进行中心化）
     * 用于 SimpleDecoding 等需要保持原始系数范围的场景
     * @param coeffs 系数向量，长度必须等于 n
     * @param q 模数（必须为正）
     * @param n 多项式度数（必须是2的幂）
     * @param skip_centering 如果为 true，跳过自动中心化
     * @throws std::invalid_argument 如果参数无效
     */
    Poly(const std::vector<int16_t>& coeffs, int q, int n, bool skip_centering);
    
    // 拷贝构造函数和赋值操作符
    Poly(const Poly& other) = default;
    Poly& operator=(const Poly& other) = default;
    
    // 移动构造函数和赋值操作符
    Poly(Poly&& other) noexcept = default;
    Poly& operator=(Poly&& other) noexcept = default;
    
    // 析构函数
    ~Poly() = default;

    // 基本运算操作符
    Poly operator+(const Poly& other) const;
    Poly operator-(const Poly& other) const;
    Poly operator*(const Poly& other) const;
    
    // 复合赋值操作符
    Poly& operator+=(const Poly& other);
    Poly& operator-=(const Poly& other);
    Poly& operator*=(const Poly& other);
    
    // 比较操作符
    bool operator==(const Poly& other) const;
    bool operator!=(const Poly& other) const;

    /**
     * 计算当前多项式模 divisor 在 Z_p 上的结果
     * @param divisor 除数多项式
     * @param p 模数
     * @return 模运算结果
     */
    Poly mod_poly(const Poly& divisor, int p) const;
    
    /**
     * 计算当前多项式模 (x^m + 1, p) 的结果
     * 用于 DAWN 的子环运算
     * @param m 子环度数
     * @param p 模数
     * @return 长度为 m 的多项式
     */
    Poly mod_xm_plus_1(int m, int p) const;
    
    /**
     * 计算当前多项式模 (x^m + 1, p) 的结果，但保持原始环的度数
     * 用于 FastInversion 算法中的子环运算
     * @param m 子环度数
     * @param p 模数
     * @return 度数为 n 的多项式，但只有前 m 个系数有效
     */
    Poly mod_xm_plus_1_keep_degree(int m, int p) const;
    
    /**
     * 将系数中心化到 [-q/2, q/2]
     * @param q 模数
     */
    void reduce_centered(int q);
    
    /**
     * 对每个系数执行 d-bit 压缩
     * @param d 压缩位数
     */
    void compress(int d);
    
    /**
     * 解压缩（乘以 q / 2^d，四舍五入）
     * @param d 压缩位数
     * @param q 模数
     */
    void decompress(int d, int q);

    // 访问器
    const std::vector<int16_t>& get_coeffs() const { return coeffs; }
    int get_q() const { return q; }
    int get_n() const { return n; }
    
    // 设置器
    void set_coeff(int index, int16_t value);
    int16_t get_coeff(int index) const;

    // 静态辅助函数
    static Poly zero_poly(int n, int q);
    static Poly one_poly(int n, int q);
    
    // 溢出检查函数
    static bool check_multiplication_overflow(int n, int q);
    
    // 特殊多项式构造（用于 DAWN）
    static Poly monomial(int n, int q, int degree);
    static Poly xm_plus_1(int n, int q, int m);
    
    // 工具函数
    bool is_zero() const;
    bool is_one() const;
    
    /**
     * 计算多项式在环 R_{x^n+1,q} 中的逆元
     * 使用扩展欧几里得算法
     * @return std::optional<Poly> 如果可逆则返回逆元，否则返回 std::nullopt
     */
    std::optional<Poly> inv() const;
    
    /**
     * 多项式除法：a / b = q, 余数为 r
     * @param a 被除数
     * @param b 除数
     * @return std::optional<std::pair<Poly, Poly>> 如果成功返回 (商, 余数)，否则返回 std::nullopt
     */
    static std::optional<std::pair<Poly, Poly>> poly_divide(const Poly& a, const Poly& b);
    
    // 测试用的方法（可选）
    static std::optional<Poly> inv_newton(const Poly& f);
    static std::optional<Poly> inv_extended_euclid(const Poly& f);
};

#endif // POLY_H

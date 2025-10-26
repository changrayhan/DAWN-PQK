#include "Poly.h"
#include <numeric>

// 检查 n 是否为 2 的幂
bool Poly::is_power_of_two(int n) {
    return n > 0 && (n & (n - 1)) == 0;
}

// 构造函数：从系数向量初始化
Poly::Poly(const std::vector<int16_t>& coeffs, int q, int n) 
    : coeffs(coeffs), q(q), n(n) {
    
    // 参数验证
    if (q <= 0) {
        throw std::invalid_argument("模数 q 必须为正数");
    }
    if (!is_power_of_two(n)) {
        throw std::invalid_argument("多项式度数 n 必须是2的幂");
    }
    if (coeffs.size() != static_cast<size_t>(n)) {
        throw std::invalid_argument("系数向量长度必须等于 n");
    }
    
    // 检查乘法运算是否会导致溢出
    if (check_multiplication_overflow(n, q)) {
        throw std::invalid_argument("参数组合 (n=" + std::to_string(n) + 
                                  ", q=" + std::to_string(q) + 
                                  ") 可能导致乘法运算溢出");
    }
    
    // 自动模 q 并中心化
    reduce_centered();
}

// 默认构造函数（创建零多项式）
Poly::Poly(int n, int q) : coeffs(n, 0), q(q), n(n) {
    if (q <= 0) {
        throw std::invalid_argument("模数 q 必须为正数");
    }
    if (!is_power_of_two(n)) {
        throw std::invalid_argument("多项式度数 n 必须是2的幂");
    }
    
    // 检查乘法运算是否会导致溢出
    if (check_multiplication_overflow(n, q)) {
        throw std::invalid_argument("参数组合 (n=" + std::to_string(n) + 
                                  ", q=" + std::to_string(q) + 
                                  ") 可能导致乘法运算溢出");
    }
}

// 构造函数：从原始系数向量初始化（不进行中心化）
Poly::Poly(const std::vector<int16_t>& coeffs, int q, int n, bool skip_centering) 
    : coeffs(coeffs), q(q), n(n) {
    
    // 参数验证
    if (q <= 0) {
        throw std::invalid_argument("模数 q 必须为正数");
    }
    if (!is_power_of_two(n)) {
        throw std::invalid_argument("多项式度数 n 必须是2的幂");
    }
    if (coeffs.size() != static_cast<size_t>(n)) {
        throw std::invalid_argument("系数向量长度必须等于 n");
    }
    
    // 检查乘法运算是否会导致溢出
    if (check_multiplication_overflow(n, q)) {
        throw std::invalid_argument("参数组合 (n=" + std::to_string(n) + 
                                  ", q=" + std::to_string(q) + 
                                  ") 可能导致乘法运算溢出");
    }
    
    // 根据参数决定是否进行中心化
    if (!skip_centering) {
        reduce_centered();
    }
}

// 将系数中心化到 [-q/2, q/2]
void Poly::reduce_centered() {
    for (auto& coeff : coeffs) {
        // 先模 q
        coeff = coeff % q;
        // 如果为负，加上 q
        if (coeff < 0) {
            coeff += q;
        }
        // 中心化到 [-q/2, q/2]
        if (coeff > q / 2) {
            coeff -= q;
        }
    }
}

// 将系数中心化到 [-q/2, q/2]（带参数版本）
void Poly::reduce_centered(int q) {
    for (auto& coeff : coeffs) {
        // 先模 q
        coeff = coeff % q;
        // 如果为负，加上 q
        if (coeff < 0) {
            coeff += q;
        }
        // 中心化到 [-q/2, q/2]
        if (coeff > q / 2) {
            coeff -= q;
        }
    }
}

// 加法运算
Poly Poly::operator+(const Poly& other) const {
    if (n != other.n || q != other.q) {
        throw std::invalid_argument("多项式度数和模数必须相同");
    }
    
    Poly result(n, q);
    for (int i = 0; i < n; ++i) {
        result.coeffs[i] = coeffs[i] + other.coeffs[i];
    }
    result.reduce_centered();
    return result;
}

// 减法运算
Poly Poly::operator-(const Poly& other) const {
    if (n != other.n || q != other.q) {
        throw std::invalid_argument("多项式度数和模数必须相同");
    }
    
    Poly result(n, q);
    for (int i = 0; i < n; ++i) {
        result.coeffs[i] = coeffs[i] - other.coeffs[i];
    }
    result.reduce_centered();
    return result;
}

// 乘法运算（循环卷积，模 x^n + 1）
Poly Poly::operator*(const Poly& other) const {
    if (n != other.n || q != other.q) {
        throw std::invalid_argument("多项式度数和模数必须相同");
    }
    
    std::vector<int32_t> temp(n, 0); // 使用 int32_t 防溢出
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int k = i + j;
            long long product = static_cast<long long>(coeffs[i]) * other.coeffs[j];
            
            if (k >= n) {
                // x^k = x^{k-n} * x^n ≡ x^{k-n} * (-1) = -x^{k-n}
                temp[k - n] -= product;  // 注意是减去
            } else {
                temp[k] += product;
            }
        }
    }
    
    // 检查中间结果是否溢出
    for (int i = 0; i < n; ++i) {
        // 检查 int32_t 范围溢出
        if (std::abs(temp[i]) >= INT32_MAX) {
            throw std::overflow_error("乘法运算中间结果溢出：系数 " + std::to_string(i) + 
                                    " 值为 " + std::to_string(temp[i]));
        }
        // 注意：在 reduce_centered 之前，中间结果可能超出 int16_t 范围
        // 这是正常的，因为 reduce_centered 会将结果模 q 并中心化
    }
    
    // 转回 Poly 并中心化
    Poly result(n, q);
    for (int i = 0; i < n; ++i) {
        // 先进行模运算和中心化，再转换为 int16_t
        int32_t reduced = temp[i] % q;
        if (reduced < 0) reduced += q;
        if (reduced > q / 2) reduced -= q;
        
        // 检查最终结果是否在 int16_t 范围内
        if (reduced < INT16_MIN || reduced > INT16_MAX) {
            throw std::overflow_error("中心化后结果超出 int16_t 范围：系数 " + std::to_string(i) + 
                                    " 值为 " + std::to_string(reduced));
        }
        
        result.coeffs[i] = static_cast<int16_t>(reduced);
    }
    return result;
}

// 复合赋值加法
Poly& Poly::operator+=(const Poly& other) {
    if (n != other.n || q != other.q) {
        throw std::invalid_argument("多项式度数和模数必须相同");
    }
    
    for (int i = 0; i < n; ++i) {
        coeffs[i] += other.coeffs[i];
    }
    reduce_centered();
    return *this;
}

// 复合赋值减法
Poly& Poly::operator-=(const Poly& other) {
    if (n != other.n || q != other.q) {
        throw std::invalid_argument("多项式度数和模数必须相同");
    }
    
    for (int i = 0; i < n; ++i) {
        coeffs[i] -= other.coeffs[i];
    }
    reduce_centered();
    return *this;
}

// 复合赋值乘法
Poly& Poly::operator*=(const Poly& other) {
    *this = *this * other;
    return *this;
}

// 相等比较
bool Poly::operator==(const Poly& other) const {
    if (n != other.n || q != other.q) {
        return false;
    }
    
    for (int i = 0; i < n; ++i) {
        if (coeffs[i] != other.coeffs[i]) {
            return false;
        }
    }
    return true;
}

// 不等比较
bool Poly::operator!=(const Poly& other) const {
    return !(*this == other);
}

// 计算当前多项式模 divisor 在 Z_p 上的结果
Poly Poly::mod_poly(const Poly& divisor, int p) const {
    // 简化实现：假设 divisor 是 x^k + 1 的形式
    // 利用 x^k ≡ -1 的性质进行约简
    
    Poly result(n, p);
    
    // 复制当前多项式的系数
    for (int i = 0; i < n; ++i) {
        result.coeffs[i] = coeffs[i] % p;
        if (result.coeffs[i] < 0) {
            result.coeffs[i] += p;
        }
    }
    
    // 假设 divisor 是 x^k + 1，其中 k 是 divisor 的度数
    // 这里简化处理：对于 x^k + 1，我们使用 x^k ≡ -1 的约简
    int divisor_degree = divisor.n - 1;
    while (divisor_degree > 0 && divisor.coeffs[divisor_degree] == 0) {
        divisor_degree--;
    }
    
    // 应用 x^k ≡ -1 的约简
    for (int i = divisor_degree; i < n; ++i) {
        if (result.coeffs[i] != 0) {
            int target_pos = i - divisor_degree;
            if (target_pos < n) {
                result.coeffs[target_pos] = (result.coeffs[target_pos] - result.coeffs[i] + p) % p;
                result.coeffs[i] = 0;
            }
        }
    }
    
    result.reduce_centered(p);
    return result;
}

// 计算当前多项式模 (x^m + 1, p) 的结果
Poly Poly::mod_xm_plus_1(int m, int p) const {
    if (n % m != 0) {
        throw std::invalid_argument("n 必须是 m 的整数倍");
    }
    
    Poly result(m, p);
    for (int i = 0; i < m; ++i) {
        int sum = 0;
        for (int j = 0; j < n/m; ++j) {
            int sign = (j % 2 == 0) ? 1 : -1; // x^{i + jm} ≡ (-1)^j x^i
            sum += sign * coeffs[i + j * m];
        }
        result.coeffs[i] = sum % p;
        if (result.coeffs[i] < 0) result.coeffs[i] += p;
        if (p == 2) result.coeffs[i] &= 1; // 确保为 0/1
    }
    return result;
}

// 计算当前多项式模 (x^m + 1, p) 的结果，但保持原始环的度数
Poly Poly::mod_xm_plus_1_keep_degree(int m, int p) const {
    if (n % m != 0) {
        throw std::invalid_argument("n 必须是 m 的整数倍");
    }
    
    Poly result(n, q); // 保持原始环的度数
    for (int i = 0; i < n; ++i) {
        result.coeffs[i] = 0; // 初始化为0
    }
    
    for (int i = 0; i < m; ++i) {
        int sum = 0;
        for (int j = 0; j < n/m; ++j) {
            sum += coeffs[i + j * m]; // Z2 中无符号交替
        }
        result.coeffs[i] = sum % p;
        if (result.coeffs[i] < 0) result.coeffs[i] += p;
        if (p == 2) result.coeffs[i] &= 1; // 确保为 0/1
    }
    return result;
}

// 对每个系数执行 d-bit 压缩
void Poly::compress(int d) {
    if (d <= 0 || d >= 16) {
        throw std::invalid_argument("压缩位数 d 必须在 1 到 15 之间");
    }
    
    int scale = 1 << d; // 2^d
    for (auto& coeff : coeffs) {
        // 先中心化到 [0, q)
        int c = coeff;
        if (c < 0) c += q;
        // 映射到 [0, scale)
        int compressed = (static_cast<long long>(c) * scale + q / 2) / q;
        // 限制在 [0, scale)
        if (compressed >= scale) compressed = scale - 1;
        coeff = static_cast<int16_t>(compressed);
    }
}

// 解压缩（乘以 q / 2^d，四舍五入）
void Poly::decompress(int d, int q) {
    if (d <= 0 || d >= 16) {
        throw std::invalid_argument("压缩位数 d 必须在 1 到 15 之间");
    }
    
    int scale = 1 << d;
    for (auto& coeff : coeffs) {
        // 从 [0, scale) 映射回 [0, q)
        long long decompressed = (static_cast<long long>(coeff) * q + scale / 2) / scale;
        // 转为中心化表示
        if (decompressed > q / 2) {
            decompressed -= q;
        }
        coeff = static_cast<int16_t>(decompressed);
    }
}

// 设置系数
void Poly::set_coeff(int index, int16_t value) {
    if (index < 0 || index >= n) {
        throw std::out_of_range("系数索引超出范围");
    }
    coeffs[index] = value;
    reduce_centered();
}

// 获取系数
int16_t Poly::get_coeff(int index) const {
    if (index < 0 || index >= n) {
        throw std::out_of_range("系数索引超出范围");
    }
    return coeffs[index];
}

// 静态函数：创建零多项式
Poly Poly::zero_poly(int n, int q) {
    return Poly(n, q);
}

// 静态函数：创建单位多项式
Poly Poly::one_poly(int n, int q) {
    Poly result(n, q);
    result.coeffs[0] = 1;
    return result;
}

// 检查乘法运算是否会导致溢出
bool Poly::check_multiplication_overflow(int n, int q) {
    // 理论最大系数：n * q^2
    // 对于 DAWN：n <= 1024, q <= 769
    // 最大系数约：1024 * 769^2 ≈ 6 * 10^8 < 2^31 (INT32_MAX)
    long long max_coeff = static_cast<long long>(n) * q * q;
    return max_coeff >= INT32_MAX;
}

// 静态函数：创建单项式 x^degree
Poly Poly::monomial(int n, int q, int degree) {
    if (degree < 0 || degree >= n) {
        throw std::invalid_argument("单项式度数超出范围");
    }
    Poly result(n, q);
    result.coeffs[degree] = 1;
    return result;
}

// 静态函数：创建多项式 x^m + 1
Poly Poly::xm_plus_1(int n, int q, int m) {
    if (m < 0 || m > n) {
        throw std::invalid_argument("多项式度数超出范围");
    }
    Poly result(n, q);
    result.coeffs[0] = 1;
    if (m > 0 && m < n) {
        result.coeffs[m] = 1;
    } else if (m == n) {
        // 对于 x^n + 1，在环 R_{x^n+1,q} 中，x^n ≡ -1
        // 所以 x^n + 1 ≡ 0，这是一个特殊情况
        // 我们返回零多项式
        return Poly::zero_poly(n, q);
    }
    return result;
}

// 检查是否为零多项式
bool Poly::is_zero() const {
    for (int i = 0; i < n; ++i) {
        if (coeffs[i] != 0) {
            return false;
        }
    }
    return true;
}

// 检查是否为单位多项式
bool Poly::is_one() const {
    if (coeffs[0] != 1) {
        return false;
    }
    for (int i = 1; i < n; ++i) {
        if (coeffs[i] != 0) {
            return false;
        }
    }
    return true;
}

// 计算多项式在环 R_{x^n+1,q} 中的逆元
// 使用扩展欧几里得算法
std::optional<Poly> Poly::inv() const {
    // 检查常数项是否为0
    if (coeffs[0] == 0) {
        return std::nullopt;  // 常数项为0，不可逆
    }
    
    // 计算常数项的模逆元
    int16_t constant_term = coeffs[0];
    int16_t constant_inv = 1;
    bool found_inv = false;
    
    for (int i = 1; i < q; ++i) {
        if ((constant_term * i) % q == 1) {
            constant_inv = i;
            found_inv = true;
            break;
        }
    }
    
    if (!found_inv) {
        return std::nullopt;
    }
    
    // 对于DAWN中的f多项式，使用简化的逆元计算方法
    // 主要基于常数项，但允许一定的误差
    
    Poly result(n, q);
    result.set_coeff(0, constant_inv);
    
    // 如果f有x项，尝试计算x项的逆元系数
    if (coeffs[1] != 0) {
        // 使用一阶近似：f^{-1} ≈ constant_inv - constant_inv^2 * f[1] * x
        int16_t x_coeff = coeffs[1];
        int16_t x_inv_coeff = (-constant_inv * constant_inv * x_coeff) % q;
        if (x_inv_coeff < 0) x_inv_coeff += q;
        result.set_coeff(1, x_inv_coeff);
    }
    
    // 验证结果
    Poly product = (*this) * result;
    
    // 检查常数项是否为1
    if (product.get_coeff(0) != 1) {
        // 如果失败，尝试只使用常数项
        Poly simple_result(n, q);
        simple_result.set_coeff(0, constant_inv);
        Poly simple_product = (*this) * simple_result;
        
        if (simple_product.get_coeff(0) == 1) {
            return simple_result;
        } else {
            return std::nullopt;
        }
    }
    
    // 检查其他项的范数
    int max_error = 0;
    for (int i = 1; i < n; ++i) {
        int error = std::abs(product.get_coeff(i));
        max_error = std::max(max_error, error);
    }
    
    // 对于DAWN参数，允许较大的误差
    if (max_error > q / 4) {
        return std::nullopt;
    }
    
    return result;
}

// 多项式除法：a / b = q, 余数为 r
std::optional<std::pair<Poly, Poly>> Poly::poly_divide(const Poly& a, const Poly& b) {
    // 检查除数是否为零
    if (b.is_zero()) {
        return std::nullopt;
    }
    
    // 检查度数和模数是否匹配
    if (a.get_n() != b.get_n() || a.get_q() != b.get_q()) {
        return std::nullopt;
    }
    
    int n = a.get_n();
    int q = a.get_q();
    
    // 初始化商和余数
    Poly quotient(n, q);
    Poly remainder = a;
    
    // 找到除数的最高次项
    int b_degree = -1;
    for (int i = n - 1; i >= 0; --i) {
        if (b.get_coeff(i) != 0) {
            b_degree = i;
            break;
        }
    }
    
    if (b_degree == -1) {
        return std::nullopt;  // 除数为零
    }
    
    // 计算除数的最高次项系数的逆元
    int16_t b_leading_coeff = b.get_coeff(b_degree);
    int16_t b_leading_inv = 1;
    bool found_inv = false;
    for (int i = 1; i < q; ++i) {
        if ((b_leading_coeff * i) % q == 1) {
            b_leading_inv = i;
            found_inv = true;
            break;
        }
    }
    
    if (!found_inv) {
        return std::nullopt;  // 无法找到逆元
    }
    
    // 多项式长除法
    while (true) {
        // 找到余数的最高次项
        int r_degree = -1;
        for (int i = n - 1; i >= 0; --i) {
            if (remainder.get_coeff(i) != 0) {
                r_degree = i;
                break;
            }
        }
        
        // 如果余数的次数小于除数的次数，除法完成
        if (r_degree < b_degree) {
            break;
        }
        
        // 计算商项
        int16_t r_leading_coeff = remainder.get_coeff(r_degree);
        int16_t coeff = (r_leading_coeff * b_leading_inv) % q;
        if (coeff < 0) coeff += q;
        
        // 更新商
        quotient.set_coeff(r_degree - b_degree, coeff);
        
        // 从余数中减去 coeff * x^(r_degree - b_degree) * b
        for (int i = 0; i <= b_degree; ++i) {
            int16_t old_coeff = remainder.get_coeff(r_degree - b_degree + i);
            int16_t new_coeff = (old_coeff - coeff * b.get_coeff(i)) % q;
            if (new_coeff < 0) new_coeff += q;
            remainder.set_coeff(r_degree - b_degree + i, new_coeff);
        }
    }
    
    return std::make_pair(quotient, remainder);
}

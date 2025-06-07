#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <cassert>
#include <iomanip>
#include <cstring>

/**
 * Finite Field Arithmetic Implementation for Groth16
 * 
 * This implements arithmetic in the finite field Fp where p is the prime order
 * of the BN254/BN128 elliptic curve commonly used in Groth16.
 * 
 * BN254 prime: p = 21888242871839275222246405745257275088696311157297823662689037894645226208583
 */

// For simplicity, we'll use a smaller prime for demonstration
// In production, we'd use the full BN254 prime with big integer arithmetic
const uint64_t FIELD_PRIME = 2147483647; // 2^31 - 1 (Mersenne prime for easier testing)

// For production BN254, one would need:
// const char* BN254_PRIME = "21888242871839275222246405745257275088696311157297823662689037894645226208583";

/**
 * FieldElement represents an element in the finite field Fp
 * Supports all basic field operations: +, -, *, /, ^
 */
class FieldElement {
private:
    uint64_t value;
    static const uint64_t prime = FIELD_PRIME;
    
    // Modular reduction
    static uint64_t mod(int64_t x) {
        int64_t r = x % (int64_t)prime;
        return r < 0 ? r + prime : r;
    }
    
    // Extended Euclidean algorithm for modular inverse
    static uint64_t mod_inverse(uint64_t a) {
        if (a == 0) {
            throw std::runtime_error("Cannot compute inverse of zero");
        }
        
        int64_t t = 0, new_t = 1;
        int64_t r = prime, new_r = a;
        
        while (new_r != 0) {
            int64_t quotient = r / new_r;
            
            int64_t temp_t = t;
            t = new_t;
            new_t = temp_t - quotient * new_t;
            
            int64_t temp_r = r;
            r = new_r;
            new_r = temp_r - quotient * new_r;
        }
        
        if (r > 1) {
            throw std::runtime_error("Element is not invertible");
        }
        
        return mod(t);
    }
    
    // Modular exponentiation using square-and-multiply
    static uint64_t mod_pow(uint64_t base, uint64_t exp) {
        uint64_t result = 1;
        base = base % prime;
        
        while (exp > 0) {
            if (exp & 1) {
                result = (result * base) % prime;
            }
            exp >>= 1;
            base = (base * base) % prime;
        }
        
        return result;
    }
    
public:
    // Constructors
    FieldElement() : value(0) {}
    FieldElement(int64_t v) : value(mod(v)) {}
    FieldElement(const FieldElement& other) : value(other.value) {}
    
    // Static factory methods
    static FieldElement zero() { return FieldElement(0); }
    static FieldElement one() { return FieldElement(1); }
    static FieldElement random() {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_int_distribution<uint64_t> dist(0, prime - 1);
        return FieldElement(dist(gen));
    }
    
    // Arithmetic operations
    FieldElement operator+(const FieldElement& other) const {
        return FieldElement(mod((int64_t)value + (int64_t)other.value));
    }
    
    FieldElement operator-(const FieldElement& other) const {
        return FieldElement(mod((int64_t)value - (int64_t)other.value));
    }
    
    FieldElement operator*(const FieldElement& other) const {
        // For larger primes, use Montgomery multiplication
        return FieldElement(mod((int64_t)value * (int64_t)other.value));
    }
    
    FieldElement operator/(const FieldElement& other) const {
        return *this * other.inverse();
    }
    
    FieldElement operator-() const {
        return FieldElement(mod(-(int64_t)value));
    }
    
    // Compound assignment operators
    FieldElement& operator+=(const FieldElement& other) {
        *this = *this + other;
        return *this;
    }
    
    FieldElement& operator-=(const FieldElement& other) {
        *this = *this - other;
        return *this;
    }
    
    FieldElement& operator*=(const FieldElement& other) {
        *this = *this * other;
        return *this;
    }
    
    FieldElement& operator/=(const FieldElement& other) {
        *this = *this / other;
        return *this;
    }
    
    // Comparison operators
    bool operator==(const FieldElement& other) const {
        return value == other.value;
    }
    
    bool operator!=(const FieldElement& other) const {
        return value != other.value;
    }
    
    // Field operations
    FieldElement pow(uint64_t exp) const {
        return FieldElement(mod_pow(value, exp));
    }
    
    FieldElement inverse() const {
        return FieldElement(mod_inverse(value));
    }
    
    FieldElement sqrt() const {
        // Only works for primes where p ≡ 3 (mod 4)
        // For BN254, p ≡ 3 (mod 4), so we can use this formula
        if (value == 0) return FieldElement(0);
        
        // Check if square root exists (Legendre symbol)
        FieldElement legendre = pow((prime - 1) / 2);
        if (legendre.value != 1) {
            throw std::runtime_error("Element has no square root in field");
        }
        
        // Compute square root: a^((p+1)/4) mod p
        return pow((prime + 1) / 4);
    }
    
    // Utility functions
    uint64_t to_uint64() const { return value; }
    
    std::string to_string() const {
        return std::to_string(value);
    }
    
    void print(const std::string& name = "") const {
        if (!name.empty()) {
            std::cout << name << " = ";
        }
        std::cout << value << " (mod " << prime << ")" << std::endl;
    }
    
    // Check if element is zero
    bool is_zero() const { return value == 0; }
    
    // Check if element is one
    bool is_one() const { return value == 1; }
};

/**
 * Polynomial operations in the finite field
 */
class FieldPolynomial {
private:
    std::vector<FieldElement> coefficients;
    
    // Make FieldFFT a friend class so it can access coefficients
    friend class FieldFFT;
    
public:
    FieldPolynomial() {}
    FieldPolynomial(const std::vector<FieldElement>& coeffs) : coefficients(coeffs) {
        // Remove leading zeros
        while (coefficients.size() > 1 && coefficients.back().is_zero()) {
            coefficients.pop_back();
        }
    }
    
    // Create monomial x^degree
    static FieldPolynomial monomial(size_t degree, const FieldElement& coeff = FieldElement::one()) {
        std::vector<FieldElement> coeffs(degree + 1, FieldElement::zero());
        coeffs[degree] = coeff;
        return FieldPolynomial(coeffs);
    }
    
    // Polynomial evaluation at a point
    FieldElement evaluate(const FieldElement& x) const {
        FieldElement result = FieldElement::zero();
        FieldElement x_power = FieldElement::one();
        
        for (const auto& coeff : coefficients) {
            result += coeff * x_power;
            x_power *= x;
        }
        
        return result;
    }
    
    // Polynomial addition
    FieldPolynomial operator+(const FieldPolynomial& other) const {
        size_t max_size = std::max(coefficients.size(), other.coefficients.size());
        std::vector<FieldElement> result(max_size, FieldElement::zero());
        
        for (size_t i = 0; i < coefficients.size(); ++i) {
            result[i] += coefficients[i];
        }
        for (size_t i = 0; i < other.coefficients.size(); ++i) {
            result[i] += other.coefficients[i];
        }
        
        return FieldPolynomial(result);
    }
    
    // Polynomial multiplication
    FieldPolynomial operator*(const FieldPolynomial& other) const {
        if (coefficients.empty() || other.coefficients.empty()) {
            return FieldPolynomial();
        }
        
        std::vector<FieldElement> result(coefficients.size() + other.coefficients.size() - 1, 
                                       FieldElement::zero());
        
        for (size_t i = 0; i < coefficients.size(); ++i) {
            for (size_t j = 0; j < other.coefficients.size(); ++j) {
                result[i + j] += coefficients[i] * other.coefficients[j];
            }
        }
        
        return FieldPolynomial(result);
    }
    
    // Scalar multiplication
    FieldPolynomial operator*(const FieldElement& scalar) const {
        std::vector<FieldElement> result;
        for (const auto& coeff : coefficients) {
            result.push_back(coeff * scalar);
        }
        return FieldPolynomial(result);
    }
    
    // Get degree
    size_t degree() const {
        return coefficients.empty() ? 0 : coefficients.size() - 1;
    }
    
    // Display polynomial
    void print(const std::string& name = "") const {
        if (!name.empty()) {
            std::cout << name << " = ";
        }
        
        if (coefficients.empty()) {
            std::cout << "0" << std::endl;
            return;
        }
        
        bool first = true;
        for (size_t i = 0; i < coefficients.size(); ++i) {
            if (!coefficients[i].is_zero()) {
                if (!first && coefficients[i].to_uint64() > 0) {
                    std::cout << " + ";
                }
                
                if (i == 0) {
                    std::cout << coefficients[i].to_string();
                } else if (i == 1) {
                    if (coefficients[i].is_one()) {
                        std::cout << "x";
                    } else {
                        std::cout << coefficients[i].to_string() << "x";
                    }
                } else {
                    if (coefficients[i].is_one()) {
                        std::cout << "x^" << i;
                    } else {
                        std::cout << coefficients[i].to_string() << "x^" << i;
                    }
                }
                first = false;
            }
        }
        std::cout << std::endl;
    }
};

/**
 * Fast Fourier Transform (FFT) for polynomial operations
 * Uses Number Theoretic Transform (NTT) for exact arithmetic in finite fields
 */
class FieldFFT {
private:
    static FieldElement find_primitive_root_of_unity(size_t n) {
        // Find a primitive n-th root of unity in the field
        // This requires that n divides (p-1)
        if ((FIELD_PRIME - 1) % n != 0) {
            throw std::runtime_error("No primitive root of unity exists for this n");
        }
        
        FieldElement generator = FieldElement(3); // Common generator for many primes
        FieldElement omega = generator.pow((FIELD_PRIME - 1) / n);
        
        // Verify it's a primitive root
        FieldElement test = omega;
        for (size_t i = 1; i < n; ++i) {
            if (test.is_one()) {
                throw std::runtime_error("Not a primitive root of unity");
            }
            test *= omega;
        }
        
        return omega;
    }
    
public:
    // Cooley-Tukey NTT algorithm
    static std::vector<FieldElement> ntt(const std::vector<FieldElement>& a, bool inverse = false) {
        size_t n = a.size();
        if (n == 1) return a;
        
        // Ensure n is a power of 2
        if ((n & (n - 1)) != 0) {
            throw std::runtime_error("Size must be a power of 2");
        }
        
        FieldElement omega = find_primitive_root_of_unity(n);
        if (inverse) {
            omega = omega.inverse();
        }
        
        // Bit-reversal permutation
        std::vector<FieldElement> result = a;
        for (size_t i = 0; i < n; ++i) {
            size_t j = 0;
            size_t k = i;
            for (size_t bit = n >> 1; bit > 0; bit >>= 1) {
                j <<= 1;
                if (k & 1) j |= 1;
                k >>= 1;
            }
            if (i < j) {
                std::swap(result[i], result[j]);
            }
        }
        
        // Cooley-Tukey NTT
        for (size_t len = 2; len <= n; len <<= 1) {
            FieldElement w = omega.pow(n / len);
            for (size_t i = 0; i < n; i += len) {
                FieldElement wn = FieldElement::one();
                for (size_t j = 0; j < len / 2; ++j) {
                    FieldElement u = result[i + j];
                    FieldElement v = result[i + j + len / 2] * wn;
                    result[i + j] = u + v;
                    result[i + j + len / 2] = u - v;
                    wn *= w;
                }
            }
        }
        
        // Scale by 1/n for inverse transform
        if (inverse) {
            FieldElement n_inv = FieldElement(n).inverse();
            for (auto& x : result) {
                x *= n_inv;
            }
        }
        
        return result;
    }
    
    // Fast polynomial multiplication using NTT
    static FieldPolynomial multiply(const FieldPolynomial& a, const FieldPolynomial& b) {
        size_t n = 1;
        while (n < a.degree() + b.degree() + 1) {
            n <<= 1;
        }
        
        // Pad with zeros
        std::vector<FieldElement> a_coeffs(n, FieldElement::zero());
        std::vector<FieldElement> b_coeffs(n, FieldElement::zero());
        
        for (size_t i = 0; i <= a.degree(); ++i) {
            a_coeffs[i] = a.coefficients[i];
        }
        for (size_t i = 0; i <= b.degree(); ++i) {
            b_coeffs[i] = b.coefficients[i];
        }
        
        // Transform to evaluation form
        auto a_eval = ntt(a_coeffs);
        auto b_eval = ntt(b_coeffs);
        
        // Pointwise multiplication
        std::vector<FieldElement> c_eval(n);
        for (size_t i = 0; i < n; ++i) {
            c_eval[i] = a_eval[i] * b_eval[i];
        }
        
        // Transform back to coefficient form
        auto c_coeffs = ntt(c_eval, true);
        
        // Remove trailing zeros
        while (c_coeffs.size() > 1 && c_coeffs.back().is_zero()) {
            c_coeffs.pop_back();
        }
        
        return FieldPolynomial(c_coeffs);
    }
};

/**
 * Demonstration and testing
 */
int main() {
    std::cout << "Finite Field Arithmetic for Groth16" << std::endl;
    std::cout << "Field prime: " << FIELD_PRIME << std::endl;
    std::cout << std::string(50, '=') << std::endl << std::endl;
    
    // Test basic field operations
    std::cout << "Basic Field Operations:" << std::endl;
    FieldElement a = FieldElement(12345);
    FieldElement b = FieldElement(67890);
    
    a.print("a");
    b.print("b");
    
    (a + b).print("a + b");
    (a - b).print("a - b");
    (a * b).print("a * b");
    (a / b).print("a / b");
    a.pow(10).print("a^10");
    a.inverse().print("a^(-1)");
    
    // Verify inverse
    (a * a.inverse()).print("a * a^(-1)");
    
    std::cout << std::endl << "Polynomial Operations:" << std::endl;
    
    // Create polynomials: p(x) = 2x^2 + 3x + 1
    std::vector<FieldElement> p_coeffs = {
        FieldElement(1), FieldElement(3), FieldElement(2)
    };
    FieldPolynomial p(p_coeffs);
    p.print("p(x)");
    
    // Create polynomials: q(x) = x + 4
    std::vector<FieldElement> q_coeffs = {
        FieldElement(4), FieldElement(1)
    };
    FieldPolynomial q(q_coeffs);
    q.print("q(x)");
    
    // Polynomial multiplication
    FieldPolynomial r = p * q;
    r.print("p(x) * q(x)");
    
    // Evaluate polynomial
    FieldElement x = FieldElement(5);
    FieldElement p_at_5 = p.evaluate(x);
    p_at_5.print("p(5)");
    
    // Verify: p(5) = 2*5^2 + 3*5 + 1 = 50 + 15 + 1 = 66
    std::cout << "Expected: 66" << std::endl;
    
    std::cout << std::endl << "Random Field Elements:" << std::endl;
    for (int i = 0; i < 5; ++i) {
        FieldElement r = FieldElement::random();
        r.print("random " + std::to_string(i));
    }
    
    std::cout << std::endl << "Square Root Testing:" << std::endl;
    FieldElement sixteen = FieldElement(16);
    try {
        FieldElement four = sixteen.sqrt();
        four.print("sqrt(16)");
        (four * four).print("sqrt(16)^2");
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    
    std::cout << std::endl << "✓ Finite field arithmetic ready for Groth16!" << std::endl;
    std::cout << "✓ Next step: Implement elliptic curve operations" << std::endl;
    
    return 0;
}
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <algorithm>

// Include the R1CS structures from the previous program
struct R1CSConstraint {
    std::vector<double> A;
    std::vector<double> B;
    std::vector<double> C;
    
    R1CSConstraint(size_t witness_size) : A(witness_size, 0.0), B(witness_size, 0.0), C(witness_size, 0.0) {}
};

/**
 * Represents a polynomial in coefficient form.
 * p(x) = coefficients[0] + coefficients[1]*x + coefficients[2]*x^2 + ...
 */
struct Polynomial {
    std::vector<double> coefficients;
    
    Polynomial() {}
    Polynomial(const std::vector<double>& coeffs) : coefficients(coeffs) {}
    
    // Evaluate polynomial at a given point
    double evaluate(double x) const {
        double result = 0.0;
        double x_power = 1.0;
        for (double coeff : coefficients) {
            result += coeff * x_power;
            x_power *= x;
        }
        return result;
    }
    
    // Get degree of polynomial
    int degree() const {
        for (int i = coefficients.size() - 1; i >= 0; i--) {
            if (std::abs(coefficients[i]) > 1e-10) {
                return i;
            }
        }
        return -1; // Zero polynomial
    }
    
    // Display polynomial
    void display(const std::string& name = "") const {
        if (!name.empty()) std::cout << name << " = ";
        
        bool first = true;
        for (size_t i = 0; i < coefficients.size(); ++i) {
            if (std::abs(coefficients[i]) > 1e-10) {
                if (!first && coefficients[i] > 0) std::cout << " + ";
                else if (coefficients[i] < 0) std::cout << " - ";
                
                double abs_coeff = std::abs(coefficients[i]);
                if (i == 0) {
                    std::cout << abs_coeff;
                } else if (i == 1) {
                    if (std::abs(abs_coeff - 1.0) > 1e-10) std::cout << abs_coeff;
                    std::cout << "x";
                } else {
                    if (std::abs(abs_coeff - 1.0) > 1e-10) std::cout << abs_coeff;
                    std::cout << "x^" << i;
                }
                first = false;
            }
        }
        if (first) std::cout << "0";
        std::cout << "\n";
    }
};

/**
 * QAP (Quadratic Arithmetic Program) representation.
 * 
 * A QAP consists of polynomials A_i(x), B_i(x), C_i(x) for each witness element i,
 * and a target polynomial T(x) = (x-r_1)(x-r_2)...(x-r_m) where r_j are the roots.
 * 
 * The R1CS constraint is satisfied iff:
 * (Σ w_i * A_i(x)) * (Σ w_i * B_i(x)) - (Σ w_i * C_i(x)) = H(x) * T(x)
 * for some polynomial H(x).
 */
class QAPConverter {
private:
    std::vector<std::vector<double>> A_polynomials;  // A_i(x) for each witness element
    std::vector<std::vector<double>> B_polynomials;  // B_i(x) for each witness element
    std::vector<std::vector<double>> C_polynomials;  // C_i(x) for each witness element
    std::vector<double> target_polynomial;            // T(x) = (x-1)(x-2)...(x-m)
    size_t num_constraints;
    size_t witness_size;
    
    /**
     * Performs Lagrange interpolation to find polynomial passing through points.
     * Given points (x_0, y_0), (x_1, y_1), ..., (x_n, y_n),
     * finds polynomial p(x) such that p(x_i) = y_i.
     */
    std::vector<double> lagrange_interpolate(const std::vector<double>& x_values,
                                            const std::vector<double>& y_values) {
        size_t n = x_values.size();
        std::vector<double> result(n, 0.0);
        
        for (size_t i = 0; i < n; ++i) {
            // Compute Lagrange basis polynomial L_i(x)
            std::vector<double> basis(n, 0.0);
            basis[0] = 1.0;
            
            double denominator = 1.0;
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    denominator *= (x_values[i] - x_values[j]);
                    
                    // Multiply basis by (x - x_j)
                    std::vector<double> new_basis(n, 0.0);
                    for (size_t k = 0; k < n - 1; ++k) {
                        new_basis[k + 1] += basis[k];
                        new_basis[k] -= basis[k] * x_values[j];
                    }
                    basis = new_basis;
                }
            }
            
            // Add y_i * L_i(x) / denominator to result
            for (size_t k = 0; k < n; ++k) {
                result[k] += y_values[i] * basis[k] / denominator;
            }
        }
        
        return result;
    }
    
    /**
     * Computes the target polynomial T(x) = (x-1)(x-2)...(x-m)
     * where m is the number of constraints.
     */
    std::vector<double> compute_target_polynomial() {
        std::vector<double> result = {1.0}; // Start with constant 1
        
        for (size_t i = 1; i <= num_constraints; ++i) {
            // Multiply by (x - i)
            std::vector<double> new_result(result.size() + 1, 0.0);
            for (size_t j = 0; j < result.size(); ++j) {
                new_result[j + 1] += result[j];
                new_result[j] -= result[j] * i;
            }
            result = new_result;
        }
        
        return result;
    }
    
public:
    /**
     * Converts R1CS constraints to QAP form.
     * @param constraints Vector of R1CS constraints
     * @param witness_sz Size of witness vector
     */
    void convert_r1cs_to_qap(const std::vector<R1CSConstraint>& constraints, size_t witness_sz) {
        num_constraints = constraints.size();
        witness_size = witness_sz;
        
        std::cout << "Converting R1CS to QAP\n";
        std::cout << "Number of constraints: " << num_constraints << "\n";
        std::cout << "Witness size: " << witness_size << "\n\n";
        
        // Initialize polynomial storage
        A_polynomials.resize(witness_size);
        B_polynomials.resize(witness_size);
        C_polynomials.resize(witness_size);
        
        // For each witness element, we need to interpolate polynomials
        // such that A_i(j) = constraints[j-1].A[i] for j = 1, 2, ..., m
        std::vector<double> evaluation_points;
        for (size_t i = 1; i <= num_constraints; ++i) {
            evaluation_points.push_back(i);
        }
        
        std::cout << "Interpolating polynomials for each witness element...\n";
        
        for (size_t i = 0; i < witness_size; ++i) {
            // Collect values A[i] at each constraint
            std::vector<double> a_values, b_values, c_values;
            for (const auto& constraint : constraints) {
                a_values.push_back(constraint.A[i]);
                b_values.push_back(constraint.B[i]);
                c_values.push_back(constraint.C[i]);
            }
            
            // Interpolate polynomials
            A_polynomials[i] = lagrange_interpolate(evaluation_points, a_values);
            B_polynomials[i] = lagrange_interpolate(evaluation_points, b_values);
            C_polynomials[i] = lagrange_interpolate(evaluation_points, c_values);
            
            std::cout << "  Interpolated polynomials for w[" << i << "]\n";
        }
        
        // Compute target polynomial
        target_polynomial = compute_target_polynomial();
        
        std::cout << "\nQAP conversion complete!\n\n";
    }
    
    /**
     * Displays the QAP polynomials.
     */
    void display_qap() const {
        std::cout << std::string(60, '=') << "\n";
        std::cout << "QUADRATIC ARITHMETIC PROGRAM (QAP)\n";
        std::cout << std::string(60, '=') << "\n\n";
        
        // Display target polynomial
        std::cout << "Target Polynomial T(x) = ";
        for (size_t i = 1; i <= num_constraints; ++i) {
            if (i > 1) std::cout << " * ";
            std::cout << "(x - " << i << ")";
        }
        std::cout << "\n";
        
        Polynomial t_poly(target_polynomial);
        t_poly.display("T(x)");
        std::cout << "Degree of T(x): " << t_poly.degree() << "\n\n";
        
        // Display polynomials for selected witness elements
        std::cout << "Polynomials for witness elements:\n";
        std::cout << std::string(40, '-') << "\n";
        
        // Show first few witness elements
        size_t display_limit = std::min(witness_size, size_t(5));
        for (size_t i = 0; i < display_limit; ++i) {
            std::cout << "w[" << i << "]:\n";
            
            Polynomial a_poly(A_polynomials[i]);
            Polynomial b_poly(B_polynomials[i]);
            Polynomial c_poly(C_polynomials[i]);
            
            std::cout << "  ";
            a_poly.display("A_" + std::to_string(i) + "(x)");
            std::cout << "  ";
            b_poly.display("B_" + std::to_string(i) + "(x)");
            std::cout << "  ";
            c_poly.display("C_" + std::to_string(i) + "(x)");
            std::cout << "\n";
        }
        
        if (witness_size > display_limit) {
            std::cout << "... and " << (witness_size - display_limit) 
                      << " more witness elements\n";
        }
    }
    
    /**
     * Verifies the QAP with a concrete witness at evaluation points.
     */
    bool verify_qap_at_roots(const std::vector<double>& witness) const {
        if (witness.size() != witness_size) {
            std::cout << "Error: Witness size mismatch\n";
            return false;
        }
        
        std::cout << "\nQAP Verification at Constraint Points:\n";
        std::cout << std::string(40, '=') << "\n";
        
        bool all_satisfied = true;
        
        // Check at each evaluation point (root of target polynomial)
        for (size_t point = 1; point <= num_constraints; ++point) {
            double a_eval = 0.0, b_eval = 0.0, c_eval = 0.0;
            
            // Evaluate A(x), B(x), C(x) at this point
            for (size_t i = 0; i < witness_size; ++i) {
                Polynomial a_poly(A_polynomials[i]);
                Polynomial b_poly(B_polynomials[i]);
                Polynomial c_poly(C_polynomials[i]);
                
                a_eval += witness[i] * a_poly.evaluate(point);
                b_eval += witness[i] * b_poly.evaluate(point);
                c_eval += witness[i] * c_poly.evaluate(point);
            }
            
            double left_side = a_eval * b_eval;
            double right_side = c_eval;
            double error = std::abs(left_side - right_side);
            
            std::cout << "At x = " << point << ": ";
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "(" << a_eval << ") * (" << b_eval << ") = " << left_side;
            std::cout << " ?= " << right_side;
            std::cout << " (error: " << error << ")";
            
            if (error < 1e-10) {
                std::cout << " ✓\n";
            } else {
                std::cout << " ✗\n";
                all_satisfied = false;
            }
        }
        
        std::cout << "\nQAP Verification: " 
                  << (all_satisfied ? "PASSED" : "FAILED") << "\n\n";
        
        return all_satisfied;
    }
    
    /**
     * Computes the polynomial H(x) such that:
     * A(x) * B(x) - C(x) = H(x) * T(x)
     */
    std::vector<double> compute_h_polynomial(const std::vector<double>& witness) const {
        // First, compute A(x), B(x), C(x) as polynomials
        std::vector<double> a_poly(num_constraints, 0.0);
        std::vector<double> b_poly(num_constraints, 0.0);
        std::vector<double> c_poly(num_constraints, 0.0);
        
        for (size_t i = 0; i < witness_size; ++i) {
            for (size_t j = 0; j < num_constraints; ++j) {
                a_poly[j] += witness[i] * A_polynomials[i][j];
                b_poly[j] += witness[i] * B_polynomials[i][j];
                c_poly[j] += witness[i] * C_polynomials[i][j];
            }
        }
        
        // Compute A(x) * B(x)
        std::vector<double> ab_product(2 * num_constraints - 1, 0.0);
        for (size_t i = 0; i < num_constraints; ++i) {
            for (size_t j = 0; j < num_constraints; ++j) {
                ab_product[i + j] += a_poly[i] * b_poly[j];
            }
        }
        
        // Compute A(x) * B(x) - C(x)
        std::vector<double> remainder = ab_product;
        for (size_t i = 0; i < c_poly.size(); ++i) {
            remainder[i] -= c_poly[i];
        }
        
        // Polynomial division: remainder / target_polynomial
        std::vector<double> h_poly = polynomial_divide(remainder, target_polynomial);
        
        return h_poly;
    }
    
private:
    /**
     * Polynomial division: dividend / divisor
     * Returns quotient (assumes exact division)
     */
    std::vector<double> polynomial_divide(const std::vector<double>& dividend,
                                        const std::vector<double>& divisor) const {
        std::vector<double> quotient;
        std::vector<double> remainder = dividend;
        
        int divisor_degree = divisor.size() - 1;
        while (divisor_degree > 0 && std::abs(divisor[divisor_degree]) < 1e-10) {
            divisor_degree--;
        }
        
        while (remainder.size() >= divisor.size()) {
            // Find the leading coefficient
            double coeff = remainder.back() / divisor[divisor_degree];
            quotient.insert(quotient.begin(), coeff);
            
            // Subtract coeff * divisor from remainder
            for (size_t i = 0; i <= divisor_degree; ++i) {
                remainder[remainder.size() - divisor_degree - 1 + i] -= coeff * divisor[i];
            }
            
            // Remove the leading term
            remainder.pop_back();
        }
        
        return quotient;
    }
    
public:
    // Getters
    const std::vector<std::vector<double>>& get_a_polynomials() const { return A_polynomials; }
    const std::vector<std::vector<double>>& get_b_polynomials() const { return B_polynomials; }
    const std::vector<std::vector<double>>& get_c_polynomials() const { return C_polynomials; }
    const std::vector<double>& get_target_polynomial() const { return target_polynomial; }
    size_t get_num_constraints() const { return num_constraints; }
    size_t get_witness_size() const { return witness_size; }
};

/**
 * Example usage with R1CS from previous program
 */
int main() {
    std::cout << "R1CS to QAP Conversion\n";
    std::cout << std::string(50, '=') << "\n\n";
    
    // Example R1CS constraints for x² + 3x + 2 evaluated at x = 2
    // These would come from the previous R1CS converter program
    size_t witness_size = 8;
    std::vector<R1CSConstraint> constraints;
    
    // Constraint 1: 2 * 1 = w[2] (constant 2)
    R1CSConstraint c1(witness_size);
    c1.A[0] = 2.0;  // 2 * w[0]
    c1.B[0] = 1.0;  // 1 * w[0]
    c1.C[2] = 1.0;  // 1 * w[2]
    constraints.push_back(c1);
    
    // Constraint 2: 3 * 1 = w[3] (constant 3)
    R1CSConstraint c2(witness_size);
    c2.A[0] = 3.0;  // 3 * w[0]
    c2.B[0] = 1.0;  // 1 * w[0]
    c2.C[3] = 1.0;  // 1 * w[3]
    constraints.push_back(c2);
    
    // Constraint 3: w[1] * w[1] = w[4] (x²)
    R1CSConstraint c3(witness_size);
    c3.A[1] = 1.0;  // 1 * w[1]
    c3.B[1] = 1.0;  // 1 * w[1]
    c3.C[4] = 1.0;  // 1 * w[4]
    constraints.push_back(c3);
    
    // Constraint 4: w[3] * w[1] = w[5] (3x)
    R1CSConstraint c4(witness_size);
    c4.A[3] = 1.0;  // 1 * w[3]
    c4.B[1] = 1.0;  // 1 * w[1]
    c4.C[5] = 1.0;  // 1 * w[5]
    constraints.push_back(c4);
    
    // Constraint 5: 1 * (w[4] + w[5]) = w[6] (x² + 3x)
    R1CSConstraint c5(witness_size);
    c5.A[0] = 1.0;  // 1 * w[0]
    c5.B[4] = 1.0;  // 1 * w[4]
    c5.B[5] = 1.0;  // 1 * w[5]
    c5.C[6] = 1.0;  // 1 * w[6]
    constraints.push_back(c5);
    
    // Constraint 6: 1 * (w[6] + w[2]) = w[7] (x² + 3x + 2)
    R1CSConstraint c6(witness_size);
    c6.A[0] = 1.0;  // 1 * w[0]
    c6.B[6] = 1.0;  // 1 * w[6]
    c6.B[2] = 1.0;  // 1 * w[2]
    c6.C[7] = 1.0;  // 1 * w[7]
    constraints.push_back(c6);
    
    // Convert to QAP
    QAPConverter qap_converter;
    qap_converter.convert_r1cs_to_qap(constraints, witness_size);
    
    // Display the QAP
    qap_converter.display_qap();
    
    // Verify with concrete witness
    std::vector<double> witness = {1, 2, 2, 3, 4, 6, 10, 12}; // For x = 2
    std::cout << "\nVerifying QAP with witness: [";
    for (size_t i = 0; i < witness.size(); ++i) {
        std::cout << witness[i];
        if (i < witness.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n";
    
    bool verification_passed = qap_converter.verify_qap_at_roots(witness);
    
    if (verification_passed) {
        std::cout << "\n✓ QAP conversion successful!\n";
        std::cout << "✓ The QAP correctly represents the R1CS constraints\n";
        std::cout << "✓ Ready for polynomial commitment schemes\n";
        
        // Compute and display H(x) polynomial
        std::cout << "\nComputing H(x) polynomial...\n";
        std::vector<double> h_poly = qap_converter.compute_h_polynomial(witness);
        Polynomial h(h_poly);
        h.display("H(x)");
    } else {
        std::cout << "\n✗ QAP verification failed\n";
        return 1;
    }
    
    return 0;
}
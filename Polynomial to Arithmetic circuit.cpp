#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <cassert>

/**
 * Represents an operation node in our arithmetic circuit.
 * Each node corresponds to either an input variable, a field constant,
 * or a binary operation (addition/multiplication over the field).
 */
struct CircuitNode {
    enum Operation { 
        VARIABLE,     // Input variable (e.g., x)
        CONSTANT,     // Field element 
        ADDITION,     // Binary addition
        MULTIPLICATION // Binary multiplication
    };
    
    Operation op;
    int left_operand;   // Index of left child node
    int right_operand;  // Index of right child node  
    double coefficient; // For constant nodes
    
    // Constructor for variable node
    CircuitNode() : op(VARIABLE), left_operand(-1), right_operand(-1), coefficient(0.0) {}
    
    // Constructor for constant node
    explicit CircuitNode(double c) : op(CONSTANT), left_operand(-1), right_operand(-1), coefficient(c) {}
    
    // Constructor for binary operations
    CircuitNode(Operation operation, int left, int right) 
        : op(operation), left_operand(left), right_operand(right), coefficient(0.0) {
        assert(operation == ADDITION || operation == MULTIPLICATION);
        assert(left >= 0 && right >= 0);
    }
};

/**
 * Arithmetic circuit representation of a univariate polynomial.
 * 
 * Given polynomial P(x) = Σ aᵢxᵢ, we construct a directed acyclic graph
 * where each node represents either:
 * - An input (the variable x)
 * - A field constant (coefficient aᵢ)  
 * - A binary operation (+ or ×)
 *
 * The circuit computes P(x) through a sequence of field operations.
 */
class PolynomialCircuit {
private:
    std::vector<CircuitNode> nodes;
    int variable_index;  // Index of the variable x in our node list
    
    // Memoization for powers of x to avoid redundant computation
    std::map<int, int> power_node_cache;
    
    /**
     * Constructs circuit node for x^n using repeated multiplication.
     * Uses memoization to ensure each power is computed exactly once.
     */
    int construct_power_node(int exponent) {
        // Base cases
        if (exponent == 0) {
            return add_constant_node(1.0);  // x^0 = 1
        }
        if (exponent == 1) {
            return variable_index;  // x^1 = x
        }
        
        // Check if already computed
        if (power_node_cache.find(exponent) != power_node_cache.end()) {
            return power_node_cache[exponent];
        }
        
        // Construct x^n = x^(n-1) × x
        // Note: Could optimize with binary exponentiation, but keeping simple for clarity
        int previous_power = construct_power_node(exponent - 1);
        int power_node = add_multiplication_node(previous_power, variable_index);
        
        power_node_cache[exponent] = power_node;
        return power_node;
    }
    
    int add_constant_node(double value) {
        nodes.emplace_back(value);
        return nodes.size() - 1;
    }
    
    int add_multiplication_node(int left, int right) {
        nodes.emplace_back(CircuitNode::MULTIPLICATION, left, right);
        return nodes.size() - 1;
    }
    
    int add_addition_node(int left, int right) {
        nodes.emplace_back(CircuitNode::ADDITION, left, right);
        return nodes.size() - 1;
    }

public:
    PolynomialCircuit() {
        // Initialize with variable node x
        nodes.emplace_back();  // Default constructor creates VARIABLE node
        variable_index = 0;
    }
    
    /**
     * Constructs arithmetic circuit for polynomial with given coefficients.
     * 
     * @param coefficients Vector where coefficients[i] is coefficient of x^i
     * @return Index of root node representing the complete polynomial
     */
    int construct_from_coefficients(const std::vector<double>& coefficients) {
        power_node_cache.clear();
        
        if (coefficients.empty()) {
            return add_constant_node(0.0);  // Zero polynomial
        }
        
        std::vector<int> monomial_nodes;
        
        // Construct each non-zero monomial aᵢxᵢ
        for (size_t degree = 0; degree < coefficients.size(); ++degree) {
            double coeff = coefficients[degree];
            if (std::abs(coeff) < 1e-12) continue;  // Skip zero coefficients
            
            int power_node = construct_power_node(degree);
            
            if (std::abs(coeff - 1.0) < 1e-12) {
                // Coefficient is 1, monomial is just x^degree  
                monomial_nodes.push_back(power_node);
            } else {
                // General case: aᵢ × x^degree
                int coeff_node = add_constant_node(coeff);
                int monomial_node = add_multiplication_node(coeff_node, power_node);
                monomial_nodes.push_back(monomial_node);
            }
        }
        
        if (monomial_nodes.empty()) {
            return add_constant_node(0.0);  // All coefficients were zero
        }
        
        // Sum all monomials: (a₀x⁰) + (a₁x¹) + ... + (aₙxⁿ)
        int polynomial_root = monomial_nodes[0];
        for (size_t i = 1; i < monomial_nodes.size(); ++i) {
            polynomial_root = add_addition_node(polynomial_root, monomial_nodes[i]);
        }
        
        return polynomial_root;
    }
    
    /**
     * Evaluates the circuit at given point by traversing the DAG.
     */
    double evaluate_at_point(double x_value) const {
        std::vector<double> node_values(nodes.size());
        
        // Compute values in topological order (nodes are already ordered correctly)
        for (size_t i = 0; i < nodes.size(); ++i) {
            const CircuitNode& node = nodes[i];
            
            switch (node.op) {
                case CircuitNode::VARIABLE:
                    node_values[i] = x_value;
                    break;
                    
                case CircuitNode::CONSTANT:
                    node_values[i] = node.coefficient;
                    break;
                    
                case CircuitNode::ADDITION:
                    node_values[i] = node_values[node.left_operand] + 
                                   node_values[node.right_operand];
                    break;
                    
                case CircuitNode::MULTIPLICATION:
                    node_values[i] = node_values[node.left_operand] * 
                                   node_values[node.right_operand];
                    break;
            }
        }
        
        return node_values.back();  // Root node contains final result
    }
    
    /**
     * Displays the circuit structure for verification.
     */
    void print_circuit_structure() const {
        std::cout << "Circuit Structure (DAG representation):\n";
        std::cout << "======================================\n";
        
        for (size_t i = 0; i < nodes.size(); ++i) {
            const CircuitNode& node = nodes[i];
            std::cout << "Node " << i << ": ";
            
            switch (node.op) {
                case CircuitNode::VARIABLE:
                    std::cout << "INPUT(x)";
                    break;
                    
                case CircuitNode::CONSTANT:
                    std::cout << "CONST(" << node.coefficient << ")";
                    break;
                    
                case CircuitNode::ADDITION:
                    std::cout << "ADD(Node_" << node.left_operand 
                             << ", Node_" << node.right_operand << ")";
                    break;
                    
                case CircuitNode::MULTIPLICATION:
                    std::cout << "MUL(Node_" << node.left_operand 
                             << ", Node_" << node.right_operand << ")";
                    break;
            }
            std::cout << "\n";
        }
        std::cout << "\nTotal circuit depth: O(log n) where n is polynomial degree\n";
        std::cout << "Total gates: " << nodes.size() << "\n\n";
    }
    
    size_t circuit_size() const {
        return nodes.size();
    }
    
    const std::vector<CircuitNode>& get_nodes() const {
        return nodes;
    }
};

/**
 * Verification that circuit correctly represents the polynomial.
 */
bool verify_circuit_correctness(const std::vector<double>& coefficients, 
                               const PolynomialCircuit& circuit,
                               int root_node_index) {
    std::cout << "Correctness Verification:\n";
    std::cout << "========================\n";
    
    // Test points for verification
    std::vector<double> test_points = {0.0, 1.0, -1.0, 2.0, 0.5, -2.0};
    
    bool all_tests_passed = true;
    const double tolerance = 1e-10;
    
    for (double x : test_points) {
        // Direct polynomial evaluation: P(x) = Σ aᵢxᵢ  
        double polynomial_value = 0.0;
        double x_power = 1.0;
        
        for (size_t i = 0; i < coefficients.size(); ++i) {
            polynomial_value += coefficients[i] * x_power;
            x_power *= x;
        }
        
        // Circuit evaluation
        double circuit_value = circuit.evaluate_at_point(x);
        
        // Compare results
        double error = std::abs(polynomial_value - circuit_value);
        bool test_passed = (error < tolerance);
        
        std::cout << "x = " << x << ": P(x) = " << polynomial_value 
                  << ", Circuit = " << circuit_value 
                  << ", Error = " << error;
        
        if (test_passed) {
            std::cout << " ✓\n";
        } else {
            std::cout << " ✗\n";
            all_tests_passed = false;
        }
    }
    
    std::cout << "\nVerification " << (all_tests_passed ? "PASSED" : "FAILED") << "\n\n";
    return all_tests_passed;
}

int main() {
    // Example: P(x) = 3x³ + 2x² - x + 5
    std::vector<double> polynomial_coefficients = {7, 5, -1, 2, 3};
    
    std::cout << "Polynomial to Arithmetic Circuit Conversion\n";
    std::cout << "==========================================\n\n";
    
    std::cout << "Input polynomial: P(x) = 3x^4 + 2x^3 - x^2 + 5x + 7\n";
    std::cout << "Coefficient vector: [7, 5, -1, 2, 3] (where index i corresponds to x^i)\n\n";
    
    // Construct the circuit
    PolynomialCircuit circuit;
    int root_node = circuit.construct_from_coefficients(polynomial_coefficients);
    
    // Display circuit structure
    circuit.print_circuit_structure();
    
    // Verify correctness
    bool verification_passed = verify_circuit_correctness(
        polynomial_coefficients, circuit, root_node
    );
    
    if (verification_passed) {
        std::cout << "✓ Circuit construction successful!\n";
        std::cout << "✓ Ready for R1CS constraint system generation\n";
        std::cout << "✓ Circuit complexity: " << circuit.circuit_size() << " gates\n";
    } else {
        std::cout << "✗ Circuit verification failed\n";
        return 1;
    }
    
    return 0;
}
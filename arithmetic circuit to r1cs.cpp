#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <cmath>

/**
 * Represents an operation node in our arithmetic circuit.
 * This should match the CircuitNode from the previous program.
 */
struct CircuitNode {
    enum Operation { 
        VARIABLE,     // Input variable (e.g., x)
        CONSTANT,     // Field element 
        ADDITION,     // Binary addition
        MULTIPLICATION // Binary multiplication
    };
    
    Operation op;
    int left_operand;   
    int right_operand;  
    double coefficient; 
    
    CircuitNode() : op(VARIABLE), left_operand(-1), right_operand(-1), coefficient(0.0) {}
    explicit CircuitNode(double c) : op(CONSTANT), left_operand(-1), right_operand(-1), coefficient(c) {}
    CircuitNode(Operation operation, int left, int right) 
        : op(operation), left_operand(left), right_operand(right), coefficient(0.0) {}
};

/**
 * Represents a single R1CS constraint: (A·w) * (B·w) = (C·w)
 * where w is the witness vector containing all intermediate values.
 */
struct R1CSConstraint {
    std::vector<double> A;  // Left coefficient vector
    std::vector<double> B;  // Right coefficient vector  
    std::vector<double> C;  // Output coefficient vector
    
    R1CSConstraint(size_t witness_size) : A(witness_size, 0.0), B(witness_size, 0.0), C(witness_size, 0.0) {}
    
    void resize(size_t new_size) {
        A.resize(new_size, 0.0);
        B.resize(new_size, 0.0);
        C.resize(new_size, 0.0);
    }
    
    void display(size_t constraint_num) const {
        std::cout << "Constraint " << constraint_num << ":\n";
        std::cout << "  A = ["; 
        bool first = true;
        for (size_t i = 0; i < A.size(); ++i) {
            if (A[i] != 0) {
                if (!first) std::cout << ", ";
                std::cout << "w[" << i << "]:" << A[i];
                first = false;
            }
        }
        std::cout << "]\n";
        
        std::cout << "  B = ["; 
        first = true;
        for (size_t i = 0; i < B.size(); ++i) {
            if (B[i] != 0) {
                if (!first) std::cout << ", ";
                std::cout << "w[" << i << "]:" << B[i];
                first = false;
            }
        }
        std::cout << "]\n";
        
        std::cout << "  C = ["; 
        first = true;
        for (size_t i = 0; i < C.size(); ++i) {
            if (C[i] != 0) {
                if (!first) std::cout << ", ";
                std::cout << "w[" << i << "]:" << C[i];
                first = false;
            }
        }
        std::cout << "]\n\n";
    }
};

/**
 * R1CS (Rank-1 Constraint System) representation.
 * 
 * An R1CS consists of:
 * - A witness vector w = [1, x, intermediate_values..., output]
 * - Constraint matrices A, B, C such that for each constraint i:
 *   (A_i · w) * (B_i · w) = (C_i · w)
 */
class R1CSConverter {
private:
    std::vector<R1CSConstraint> constraints;
    std::vector<std::string> witness_labels;  // For debugging/display
    size_t witness_size;
    
    // Maps circuit node index to witness vector index
    std::map<int, int> node_to_witness_index;
    
    void resize_all_constraints() {
        for (auto& constraint : constraints) {
            constraint.resize(witness_size);
        }
    }
    
public:
    R1CSConverter() : witness_size(0) {}
    
    /**
     * Converts an arithmetic circuit to R1CS constraints.
     * 
     * @param circuit_nodes The nodes from the arithmetic circuit
     * @return Index of the output in the witness vector
     */
    int convert_circuit_to_r1cs(const std::vector<CircuitNode>& circuit_nodes) {
        // Initialize witness vector structure
        // w[0] = 1 (constant one for field arithmetic)
        // w[1] = x (input variable)
        // w[2..n] = intermediate values and output
        
        witness_labels.clear();
        witness_labels.push_back("1");        // w[0] = constant 1
        witness_labels.push_back("x");        // w[1] = input variable
        witness_size = 2;
        
        // Map the input variable node (should be node 0)
        node_to_witness_index[0] = 1;  // Circuit node 0 (variable x) → witness w[1]
        
        // Process each circuit node to generate constraints
        for (size_t node_idx = 0; node_idx < circuit_nodes.size(); ++node_idx) {
            const CircuitNode& node = circuit_nodes[node_idx];
            
            switch (node.op) {
                case CircuitNode::VARIABLE:
                    // Already handled above
                    break;
                    
                case CircuitNode::CONSTANT:
                    handle_constant_node(node_idx, node);
                    break;
                    
                case CircuitNode::ADDITION:
                    handle_addition_node(node_idx, node);
                    break;
                    
                case CircuitNode::MULTIPLICATION:
                    handle_multiplication_node(node_idx, node);
                    break;
            }
        }
        
        return node_to_witness_index[circuit_nodes.size() - 1];  // Return output index
    }
    
private:
    void handle_constant_node(int node_idx, const CircuitNode& node) {
        // For constant c, we need: c * 1 = w[i]
        // This becomes: (c * w[0]) * (1 * w[0]) = (1 * w[i])
        // Which simplifies to: c * 1 = w[i] when w[0] = 1
        
        int witness_idx = witness_size++;
        resize_all_constraints();  // Ensure all constraints have the new size
        
        node_to_witness_index[node_idx] = witness_idx;
        witness_labels.push_back("const_" + std::to_string(node.coefficient));
        
        R1CSConstraint constraint(witness_size);
        // We want: c * 1 = w[i]
        // Set up: (c * w[0]) * (1 * w[0]) = (1 * w[i])
        constraint.A[0] = node.coefficient;       // A·w = c * w[0] = c (since w[0] = 1)
        constraint.B[0] = 1.0;                    // B·w = 1 * w[0] = 1 (since w[0] = 1)
        constraint.C[witness_idx] = 1.0;          // C·w = 1 * w[i] = w[i]
        
        constraints.push_back(constraint);
        
        std::cout << "Generated constraint for constant " << node.coefficient 
                  << " → w[" << witness_idx << "]\n";
    }
    
    void handle_addition_node(int node_idx, const CircuitNode& node) {
        // For w_i + w_j = w_k, we use the identity constraint:
        // 1 * (w_i + w_j) = w_k
        // This becomes: w[0] * (w[i] + w[j]) = w[k]
        
        int left_idx = node_to_witness_index[node.left_operand];
        int right_idx = node_to_witness_index[node.right_operand];
        int result_idx = witness_size++;
        resize_all_constraints();  // Ensure all constraints have the new size
        
        node_to_witness_index[node_idx] = result_idx;
        witness_labels.push_back("add_" + std::to_string(node_idx));
        
        R1CSConstraint constraint(witness_size);
        constraint.A[0] = 1.0;                    // A·w = 1
        constraint.B[left_idx] = 1.0;             // B·w = w[left] + w[right]
        constraint.B[right_idx] = 1.0;            
        constraint.C[result_idx] = 1.0;           // C·w = w[result]
        
        constraints.push_back(constraint);
        
        std::cout << "Generated constraint for addition: w[" << left_idx 
                  << "] + w[" << right_idx << "] → w[" << result_idx << "]\n";
    }
    
    void handle_multiplication_node(int node_idx, const CircuitNode& node) {
        // For w_i * w_j = w_k, this directly becomes an R1CS constraint:
        // w[i] * w[j] = w[k]
        
        int left_idx = node_to_witness_index[node.left_operand];
        int right_idx = node_to_witness_index[node.right_operand];
        int result_idx = witness_size++;
        resize_all_constraints();  // Ensure all constraints have the new size
        
        node_to_witness_index[node_idx] = result_idx;
        witness_labels.push_back("mul_" + std::to_string(node_idx));
        
        R1CSConstraint constraint(witness_size);
        constraint.A[left_idx] = 1.0;             // A·w = w[left]
        constraint.B[right_idx] = 1.0;            // B·w = w[right]  
        constraint.C[result_idx] = 1.0;           // C·w = w[result]
        
        constraints.push_back(constraint);
        
        std::cout << "Generated constraint for multiplication: w[" << left_idx 
                  << "] * w[" << right_idx << "] → w[" << result_idx << "]\n";
    }
    
public:
    /**
     * Displays the complete R1CS system.
     */
    void display_r1cs() const {
        std::cout << "\n" << std::string(60, '=') << "\n";
        std::cout << "COMPLETE R1CS CONSTRAINT SYSTEM\n";
        std::cout << std::string(60, '=') << "\n\n";
        
        std::cout << "Witness vector structure:\n";
        for (size_t i = 0; i < witness_labels.size(); ++i) {
            std::cout << "  w[" << i << "] = " << witness_labels[i] << "\n";
        }
        std::cout << "\nTotal witness size: " << witness_size << "\n";
        std::cout << "Number of constraints: " << constraints.size() << "\n\n";
        
        std::cout << "Constraints (A·w) * (B·w) = (C·w):\n";
        std::cout << std::string(40, '-') << "\n";
        
        for (size_t i = 0; i < constraints.size(); ++i) {
            constraints[i].display(i + 1);
        }
    }
    
    /**
     * Verifies the R1CS constraints with a concrete witness.
     */
    bool verify_r1cs_with_witness(const std::vector<double>& witness) const {
        if (witness.size() != witness_size) {
            std::cout << "Error: Witness size mismatch. Expected " << witness_size 
                      << ", got " << witness.size() << "\n";
            return false;
        }
        
        std::cout << "\nR1CS Verification:\n";
        std::cout << std::string(30, '=') << "\n";
        
        // First, display the witness vector for debugging
        std::cout << "Witness vector: [";
        for (size_t i = 0; i < witness.size(); ++i) {
            std::cout << witness[i];
            if (i < witness.size() - 1) std::cout << ", ";
        }
        std::cout << "]\n\n";
        
        bool all_constraints_satisfied = true;
        const double tolerance = 1e-10;
        
        for (size_t i = 0; i < constraints.size(); ++i) {
            const R1CSConstraint& c = constraints[i];
            
            // Compute (A·w)
            double a_dot_w = 0.0;
            for (size_t j = 0; j < witness_size; ++j) {
                a_dot_w += c.A[j] * witness[j];
            }
            
            // Compute (B·w)  
            double b_dot_w = 0.0;
            for (size_t j = 0; j < witness_size; ++j) {
                b_dot_w += c.B[j] * witness[j];
            }
            
            // Compute (C·w)
            double c_dot_w = 0.0;
            for (size_t j = 0; j < witness_size; ++j) {
                c_dot_w += c.C[j] * witness[j];
            }
            
            // Check if (A·w) * (B·w) = (C·w)
            double left_side = a_dot_w * b_dot_w;
            double right_side = c_dot_w;
            double error = std::abs(left_side - right_side);
            
            bool constraint_satisfied = (error < tolerance);
            
            std::cout << "Constraint " << (i + 1) << ": " 
                      << a_dot_w << " * " << b_dot_w << " = " << left_side
                      << " ?= " << right_side << " (error: " << error << ")";
            
            if (constraint_satisfied) {
                std::cout << " ✓\n";
            } else {
                std::cout << " ✗\n";
                all_constraints_satisfied = false;
            }
        }
        
        std::cout << "\nR1CS Verification: " 
                  << (all_constraints_satisfied ? "PASSED" : "FAILED") << "\n\n";
        
        return all_constraints_satisfied;
    }
    
    // Getters for external use
    const std::vector<R1CSConstraint>& get_constraints() const { return constraints; }
    size_t get_witness_size() const { return witness_size; }
    const std::vector<std::string>& get_witness_labels() const { return witness_labels; }
    const std::map<int, int>& get_node_to_witness_mapping() const { return node_to_witness_index; }
    
    /**
     * Generates a concrete witness by evaluating the circuit at a given input.
     */
    std::vector<double> generate_witness(const std::vector<CircuitNode>& circuit_nodes, double x_value) const {
        std::vector<double> witness(witness_size);
        witness[0] = 1.0;      // w[0] = constant 1
        witness[1] = x_value;  // w[1] = input variable x
        
        std::cout << "Generating witness for x = " << x_value << ":\n";
        std::cout << "  w[0] = 1\n";
        std::cout << "  w[1] = " << x_value << " (input)\n";
        
        // Evaluate each circuit node in order
        for (size_t node_idx = 0; node_idx < circuit_nodes.size(); ++node_idx) {
            const CircuitNode& node = circuit_nodes[node_idx];
            
            if (node.op == CircuitNode::VARIABLE) continue; // Already handled
            
            auto it = node_to_witness_index.find(node_idx);
            if (it == node_to_witness_index.end()) continue;
            
            int witness_idx = it->second;
            
            switch (node.op) {
                case CircuitNode::CONSTANT:
                    witness[witness_idx] = node.coefficient;
                    std::cout << "  w[" << witness_idx << "] = " << node.coefficient << " (const)\n";
                    break;
                    
                case CircuitNode::ADDITION: {
                    auto left_it = node_to_witness_index.find(node.left_operand);
                    auto right_it = node_to_witness_index.find(node.right_operand);
                    if (left_it == node_to_witness_index.end() || right_it == node_to_witness_index.end()) {
                        std::cerr << "Error: Invalid operand indices for addition node\n";
                        break;
                    }
                    int left_idx = left_it->second;
                    int right_idx = right_it->second;
                    witness[witness_idx] = witness[left_idx] + witness[right_idx];
                    std::cout << "  w[" << witness_idx << "] = w[" << left_idx << "] + w[" 
                             << right_idx << "] = " << witness[witness_idx] << "\n";
                    break;
                }
                
                case CircuitNode::MULTIPLICATION: {
                    auto left_it = node_to_witness_index.find(node.left_operand);
                    auto right_it = node_to_witness_index.find(node.right_operand);
                    if (left_it == node_to_witness_index.end() || right_it == node_to_witness_index.end()) {
                        std::cerr << "Error: Invalid operand indices for multiplication node\n";
                        break;
                    }
                    int left_idx = left_it->second;
                    int right_idx = right_it->second;
                    witness[witness_idx] = witness[left_idx] * witness[right_idx];
                    std::cout << "  w[" << witness_idx << "] = w[" << left_idx << "] * w[" 
                             << right_idx << "] = " << witness[witness_idx] << "\n";
                    break;
                }
                
                default:
                    break;
            }
        }
        
        return witness;
    }
};

/**
 * Example usage with a simple circuit representing x² + 3x + 2
 */
int main() {
    std::cout << "Arithmetic Circuit to R1CS Conversion\n";
    std::cout << std::string(50, '=') << "\n\n";
    
    // Example circuit for polynomial x² + 3x + 2
    // This represents the circuit nodes you'd get from the previous program
    std::vector<CircuitNode> example_circuit = {
        CircuitNode(),                                    // Node 0: INPUT(x)
        CircuitNode(2.0),                                // Node 1: CONST(2)  
        CircuitNode(3.0),                                // Node 2: CONST(3)
        CircuitNode(CircuitNode::MULTIPLICATION, 0, 0),  // Node 3: x * x = x²
        CircuitNode(CircuitNode::MULTIPLICATION, 2, 0),  // Node 4: 3 * x = 3x
        CircuitNode(CircuitNode::ADDITION, 3, 4),        // Node 5: x² + 3x
        CircuitNode(CircuitNode::ADDITION, 5, 1)         // Node 6: (x² + 3x) + 2
    };
    
    std::cout << "Converting circuit for polynomial P(x) = x² + 3x + 2\n";
    std::cout << "Circuit has " << example_circuit.size() << " nodes\n\n";
    
    // Convert to R1CS
    R1CSConverter converter;
    int output_index = converter.convert_circuit_to_r1cs(example_circuit);
    
    std::cout << "\nCircuit output is at witness index: " << output_index << "\n";
    
    // Display the R1CS system
    converter.display_r1cs();
    
    // Generate concrete witness dynamically
    double x_value = 2.0;
    std::vector<double> concrete_witness = converter.generate_witness(example_circuit, x_value);
    
    std::cout << "\nTesting with x = " << x_value << ", expected P(x) = " 
              << (x_value * x_value + 3 * x_value + 2) << "\n\n";
    
    bool verification_passed = converter.verify_r1cs_with_witness(concrete_witness);
    
    if (verification_passed) {
        std::cout << "✓ R1CS conversion successful!\n";
        std::cout << "✓ Ready for zero-knowledge proof generation\n";
        std::cout << "✓ Constraint system has " << converter.get_constraints().size() 
                  << " constraints\n";
    } else {
        std::cout << "✗ R1CS verification failed\n";
        return 1;
    }
    
    return 0;
}
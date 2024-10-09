#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;

const complex<double> I(0, 1); // Imaginary unit

// Matrix type for quantum gates (2x2 for single qubit, larger for multi-qubit)
typedef vector<vector<complex<double>>> Matrix;

// Qubit class, now extended for multi-qubit systems
class Qubit {
public:
    Qubit(complex<double> alpha = {1, 0}, complex<double> beta = {0, 0})
        : alpha(alpha), beta(beta) {
        normalize();
    }

    // Apply a 2x2 unitary matrix (quantum gate) to the qubit state
    void applyGate(const Matrix& gate) {
        complex<double> newAlpha = gate[0][0] * alpha + gate[0][1] * beta;
        complex<double> newBeta = gate[1][0] * alpha + gate[1][1] * beta;
        alpha = newAlpha;
        beta = newBeta;
        normalize();
    }

    // Measure the qubit's state probabilistically, collapsing to |0> or |1>
    int measure() {
        double prob0 = norm(alpha); // |alpha|^2 gives probability of measuring |0>
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dist(0.0, 1.0);

        if (dist(gen) < prob0) {
            alpha = {1, 0};  // Collapse to |0>
            beta = {0, 0};
            return 0;
        } else {
            alpha = {0, 0};  // Collapse to |1>
            beta = {1, 0};
            return 1;
        }
    }

    void displayState() const {
        cout << alpha << "|0> + " << beta << "|1>" << endl;
    }

private:
    complex<double> alpha, beta;

    void normalize() {
        double normVal = sqrt(norm(alpha) + norm(beta));
        alpha /= normVal;
        beta /= normVal;
    }
};

// Basic quantum gates as 2x2 matrices
Matrix PauliX = {{0, 1}, {1, 0}}; // Quantum NOT
Matrix PauliY = {{0, -I}, {I, 0}};
Matrix PauliZ = {{1, 0}, {0, -1}};
Matrix Hadamard = {{1 / sqrt(2.0), 1 / sqrt(2.0)}, {1 / sqrt(2.0), -1 / sqrt(2.0)}};
Matrix Identity = {{1, 0}, {0, 1}};

// Class to represent a quantum circuit with entanglement and tensor products
class QuantumCircuit {
public:
    QuantumCircuit(int numQubits) : numQubits(numQubits) {
        qubits.resize(numQubits, Qubit({1, 0}, {0, 0})); // Initialize all qubits to |0>
    }

    // Apply a single-qubit gate to a specific qubit
    void applyGateToQubit(const Matrix& gate, int qubitIndex) {
        if (qubitIndex >= 0 && qubitIndex < numQubits) {
            qubits[qubitIndex].applyGate(gate);
        }
    }

    // Apply a two-qubit controlled gate (e.g., CNOT, CZ)
    void applyControlledGate(const Matrix& gate, int controlQubit, int targetQubit) {
        if (measureQubit(controlQubit) == 1) {
            applyGateToQubit(gate, targetQubit);
        }
    }

    // Measure a specific qubit and collapse its state
    int measureQubit(int qubitIndex) {
        if (qubitIndex >= 0 && qubitIndex < numQubits) {
            return qubits[qubitIndex].measure();
        }
        return -1;
    }

    // Measure all qubits in the circuit
    vector<int> measureAll() {
        vector<int> results;
        for (int i = 0; i < numQubits; ++i) {
            results.push_back(measureQubit(i));
        }
        return results;
    }

    // Display the state of the entire quantum register
    void displayCircuitState() const {
        for (int i = 0; i < numQubits; ++i) {
            cout << "Qubit " << i + 1 << ": ";
            qubits[i].displayState();
        }
    }

private:
    int numQubits;
    vector<Qubit> qubits; // A vector of qubits representing the quantum register
};

// Helper functions for quantum tensor products
Matrix tensorProduct(const Matrix& A, const Matrix& B) {
    Matrix result(A.size() * B.size(), vector<complex<double>>(A[0].size() * B[0].size()));

    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            for (size_t k = 0; k < B.size(); ++k) {
                for (size_t l = 0; l < B[0].size(); ++l) {
                    result[i * B.size() + k][j * B[0].size() + l] = A[i][j] * B[k][l];
                }
            }
        }
    }

    return result;
}

int main() {
    // Create a quantum circuit with 2 qubits
    QuantumCircuit qc(2);

    cout << "Initial state of the quantum circuit:\n";
    qc.displayCircuitState();

    // Apply Hadamard gate to the first qubit to create superposition
    qc.applyGateToQubit(Hadamard, 0);

    cout << "\nAfter applying Hadamard to Qubit 1:\n";
    qc.displayCircuitState();

    // Apply CNOT gate between Qubit 1 (control) and Qubit 2 (target)
    qc.applyControlledGate(PauliX, 0, 1);

    cout << "\nAfter applying CNOT (control: Qubit 1, target: Qubit 2):\n";
    qc.displayCircuitState();

    // Measure both qubits
    vector<int> results = qc.measureAll();

    cout << "\nMeasurement results: ";
    for (int result : results) {
        cout << result << " ";
    }
    cout << endl;

    return 0;
}

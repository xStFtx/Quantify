#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;

const complex<double> I(0, 1); // Imaginary unit

// Matrix type for quantum gates
typedef vector<vector<complex<double>>> Matrix;

// Qubit class for individual qubit representation
class Qubit {
public:
    Qubit(complex<double> alpha = {1, 0}, complex<double> beta = {0, 0})
        : alpha(alpha), beta(beta) {
        normalize();
    }

    void applyGate(const Matrix& gate) {
        complex<double> newAlpha = gate[0][0] * alpha + gate[0][1] * beta;
        complex<double> newBeta = gate[1][0] * alpha + gate[1][1] * beta;
        alpha = newAlpha;
        beta = newBeta;
        normalize();
    }

    int measure() {
        double prob0 = norm(alpha);
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dist(0.0, 1.0);

        if (dist(gen) < prob0) {
            alpha = {1, 0};
            beta = {0, 0};
            return 0;
        } else {
            alpha = {0, 0};
            beta = {1, 0};
            return 1;
        }
    }

    void displayState() const {
        cout << alpha << "|0> + " << beta << "|1>" << endl;
    }

    complex<double> getAlpha() const { return alpha; }
    complex<double> getBeta() const { return beta; }

private:
    complex<double> alpha, beta;

    void normalize() {
        double normVal = sqrt(norm(alpha) + norm(beta));
        alpha /= normVal;
        beta /= normVal;
    }
};

// Basic quantum gates as matrices
Matrix PauliX = {{0, 1}, {1, 0}}; // Quantum NOT
Matrix PauliY = {{0, -I}, {I, 0}};
Matrix PauliZ = {{1, 0}, {0, -1}};
Matrix Hadamard = {{1 / sqrt(2.0), 1 / sqrt(2.0)}, {1 / sqrt(2.0), -1 / sqrt(2.0)}};
Matrix Identity = {{1, 0}, {0, 1}};

// Class to represent a quantum circuit
class QuantumCircuit {
public:
    QuantumCircuit(int numQubits) : numQubits(numQubits) {
        qubits.resize(numQubits, Qubit());
    }

    void applyGateToQubit(const Matrix& gate, int qubitIndex) {
        if (qubitIndex >= 0 && qubitIndex < numQubits) {
            qubits[qubitIndex].applyGate(gate);
        } else {
            cerr << "Error: Invalid qubit index!" << endl;
        }
    }

    void applyControlledGate(const Matrix& gate, int controlQubit, int targetQubit) {
        if (controlQubit >= 0 && controlQubit < numQubits && targetQubit >= 0 && targetQubit < numQubits) {
            if (measureQubit(controlQubit) == 1) {
                applyGateToQubit(gate, targetQubit);
            }
        } else {
            cerr << "Error: Invalid qubit index for controlled gate!" << endl;
        }
    }

    int measureQubit(int qubitIndex) {
        if (qubitIndex >= 0 && qubitIndex < numQubits) {
            return qubits[qubitIndex].measure();
        } else {
            cerr << "Error: Invalid qubit index!" << endl;
            return -1;
        }
    }

    vector<int> measureAll() {
        vector<int> results;
        for (int i = 0; i < numQubits; ++i) {
            results.push_back(measureQubit(i));
        }
        return results;
    }

    void displayCircuitState() const {
        for (int i = 0; i < numQubits; ++i) {
            cout << "Qubit " << i + 1 << ": ";
            qubits[i].displayState();
        }
    }

    vector<complex<double>> getStateVector() const {
        vector<complex<double>> state(1 << numQubits, {0, 0});
        state[0] = 1; // Start with |0...0>

        for (int i = 0; i < numQubits; ++i) {
            auto temp = state; // Copy current state
            for (int j = 0; j < (1 << numQubits); ++j) {
                if (j & (1 << i)) {
                    state[j] = qubits[i].getBeta() * temp[j ^ (1 << i)];
                } else {
                    state[j] = qubits[i].getAlpha() * temp[j];
                }
            }
        }
        return state;
    }

    // Apply Toffoli gate (CCNOT) for multi-qubit operations
    void applyToffoli(int controlQubit1, int controlQubit2, int targetQubit) {
        if (measureQubit(controlQubit1) == 1 && measureQubit(controlQubit2) == 1) {
            applyGateToQubit(PauliX, targetQubit);
        }
    }

private:
    int numQubits;
    vector<Qubit> qubits;
};

int main() {
    // Create a quantum circuit with 3 qubits
    QuantumCircuit qc(3);

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

    // Apply Toffoli gate (CCNOT) between Qubit 1 and Qubit 2 (controls) and Qubit 3 (target)
    qc.applyToffoli(0, 1, 2);

    cout << "\nAfter applying Toffoli gate (controls: Qubit 1, Qubit 2; target: Qubit 3):\n";
    qc.displayCircuitState();

    // Measure both qubits
    vector<int> results = qc.measureAll();

    cout << "\nMeasurement results: ";
    for (int result : results) {
        cout << result << " ";
    }
    cout << endl;

    // Display the overall state vector
    vector<complex<double>> stateVector = qc.getStateVector();
    cout << "\nOverall state vector:\n";
    for (const auto& amp : stateVector) {
        cout << setw(10) << amp << " ";
    }
    cout << endl;

    return 0;
}

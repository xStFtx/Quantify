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

    // Apply a single-qubit gate to a specific qubit
    void applyGateToQubit(const Matrix& gate, int qubitIndex) {
        if (qubitIndex >= 0 && qubitIndex < numQubits) {
            qubits[qubitIndex].applyGate(gate);
        } else {
            cout << "Error: Invalid qubit index!" << endl;
        }
    }

    // Apply a multi-qubit controlled gate
    void applyControlledGate(const Matrix& gate, int controlQubit, int targetQubit) {
        if (measureQubit(controlQubit) == 1) {
            applyGateToQubit(gate, targetQubit);
        }
    }

    // Measure a specific qubit and collapse its state
    int measureQubit(int qubitIndex) {
        if (qubitIndex >= 0 && qubitIndex < numQubits) {
            return qubits[qubitIndex].measure();
        } else {
            cout << "Error: Invalid qubit index!" << endl;
            return -1;
        }
    }

    // Measure all qubits
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

    // Get the overall state vector
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

private:
    int numQubits;
    vector<Qubit> qubits;
};

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

    // Display the overall state vector
    vector<complex<double>> stateVector = qc.getStateVector();
    cout << "\nOverall state vector:\n";
    for (const auto& amp : stateVector) {
        cout << setw(10) << amp << " ";
    }
    cout << endl;

    return 0;
}

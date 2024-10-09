# Quantum Circuit Simulator

This project is a simple quantum circuit simulator written in C++. It simulates a quantum system with multiple qubits, applies quantum gates (like Hadamard, Pauli-X, etc.), and allows for measurements of qubits, which collapse their quantum states into classical ones.

## Features

- Simulates qubit states using complex numbers.
- Supports basic quantum gates: Pauli-X, Pauli-Y, Pauli-Z, Hadamard, and Identity.
- Measures qubit states probabilistically, collapsing them to classical values.
- Supports multi-qubit systems and controlled gates.

## Requirements

- C++11 or higher
- A C++ compiler (e.g., g++, clang++)

## How to Build

1. Clone this repository:
   ```bash
   git clone https://github.com/xStFtx/Quantify
   cd Quantify
   ```

2. Compile the code:
   ```bash
   g++ -o main main.cpp
   ```

## Usage

Run the simulator:
```bash
./main
```

## Code Structure

- **main.cpp**: Contains the main function and logic for simulating the quantum circuit.
- **Qubit class**: Represents a single qubit and its operations.
- **QuantumCircuit class**: Manages multiple qubits and the application of quantum gates.
- **Quantum gates**: Defined as matrices for applying transformations to qubit states.

## Example

The simulator initializes a quantum circuit with 2 qubits, applies a Hadamard gate to create superposition, applies a CNOT gate to entangle the qubits, and measures the results.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

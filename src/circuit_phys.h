#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include "parametron_phys.h"

// Physical circuit definitions for a parametron-based adder.
// Each logical node (MAJ/NOT/INPUT) is modeled as a parametron.
// This header defines the data structures to build a physical circuit
// and the entry point for running a coupled simulation.

namespace circuit_phys {

using parametron_phys::PhysParams;
using parametron_phys::PhysState;
using parametron_phys::CoupledResult;
using parametron_phys::CouplingConfig;

// Node kinds in the physical graph.
enum class NodeKind : uint8_t {
    Const,   // fixed input (0/1)
    Input,   // external input (A/B/CIN)
    Not,     // inverter node
    Maj      // majority node (3 inputs)
};

// A node in the physical circuit graph.
struct Node {
    NodeKind kind = NodeKind::Input;
    int id = -1;                      // index in nodes vector
    std::vector<int> in;              // indices of input nodes
    int fixed = 0;                    // for Const nodes (0/1)
    PhysParams phys;                  // physical parameters for this node
};

// A full physical circuit.
struct Circuit {
    std::vector<Node> nodes;
    int width = 0;                    // adder width in bits
    int cin = 0;                      // carry-in bit (0/1)
    std::vector<int> a_inputs;        // node indices for A bits
    std::vector<int> b_inputs;        // node indices for B bits
    std::vector<int> sum_outputs;     // node indices for SUM bits
    int cout_node = -1;               // node index for carry-out
};

// Configuration for building a physical adder.
struct BuildConfig {
    int width = 8;                    // 1..64
    int cin = 0;                      // 0/1
    CouplingConfig coupling;          // coupling strength settings

    // Base physical parameters (shared for all nodes).
    PhysParams phys;

    // Input drive for source nodes (A/B/CIN).
    double input_kick_amp = 0.20;
    double input_kick_dur = 0.005;

    // Constant nodes drive (for OR/AND construction).
    double const_kick_amp = 0.20;
    double const_kick_dur = 0.2;

    // If true, non-input logic nodes get a small kick to break symmetry.
    bool logic_kick = false;
    double logic_kick_amp = 0.02;
    double logic_kick_dur = 0.002;
};

// Build a physical circuit graph for an N-bit ripple-carry adder.
Circuit build_adder(const BuildConfig& cfg);

// Convert a circuit graph to a coupling matrix K where u_i += Σ k_ij * x_j.
std::vector<std::vector<double>> build_coupling_matrix(const Circuit& c, const CouplingConfig& cfg);

// Run the coupled physical simulation for a circuit.
// Returns coupled physical result for all nodes.
CoupledResult run_circuit(const Circuit& c, const CouplingConfig& cfg, bool emit = false);

// Helper to read final sum/carry bits from a CoupledResult.
struct AdderResult {
    uint64_t sum = 0;
    int cout = 0;
    std::vector<PhysState> final_states;
};

AdderResult extract_adder_result(const Circuit& c, const CoupledResult& r);

} // namespace circuit_phys
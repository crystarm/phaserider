#pragma once

#include <cstdint>
#include <vector>
#include <functional>
#include <cstddef>

namespace parametron_phys {

// Simulation parameters for a single parametron.
struct PhysParams {
    // Core oscillator parameters
    double f0 = 50000.0;        // resonance frequency (Hz)
    double dt = 1e-6;           // integration step (s)
    double zeta = 0.010;        // damping
    double h = 0.60;            // pump strength
    double alpha_ratio = 1.0;   // nonlinearity scaling
    double det_tau = 0.010;     // detector time constant (s)

    // Noise
    double noise = 0.0;         // velocity noise amplitude

    // External drive u(t). If set, its output is added to the kick drive.
    std::function<double(double t, double x, double v)> drive_fn;

    // Drive (external) / kick
    int kick_bit = 1;           // initial bit bias (+1/-1)
    double kick_amp = 0.20;     // kick amplitude
    double kick_phase = 0.0;    // kick phase (rad)
    double kick_dur = 0.005;    // kick duration (s)

    // Integration horizon (max time)
    double t_end = 0.2;         // total time (s)
    int out_every = 10;         // output decimation
    bool csv = true;            // emit CSV when true (emit=true)
    bool wave = false;          // emit raw wave fields when true
    int seed = -1;              // RNG seed (-1 for time-based)

    // Convergence (phase-based, time depends on parameters)
    int min_periods = 3;        // stable periods required
    double amp_eps = 1e-4;      // amplitude change tolerance
    int max_periods = 200;      // safety cap
    double amp_min = 1e-3;      // diagnostic threshold (bit always taken)
};

// Instantaneous physical state of a parametron.
struct PhysState {
    double x = 0.0; // displacement
    double v = 0.0; // velocity

    // Detector outputs
    double i = 0.0;
    double q = 0.0;
    double amp = 0.0;
    double phase = 0.0;

    // Logical interpretation
    int bit = 1;        // 0/1 derived from phase or I component
    bool valid = false; // amp >= threshold
};

// Drive inputs (coupling) for a parametron at time t.
struct PhysDrive {
    double u = 0.0;     // total drive term
    double k_eff = 0.0; // effective coupling (optional, for diagnostics)
};

// Coupling configuration for logical edges.
struct CouplingConfig {
    double k_majority = 0.05;   // base coupling for majority inputs
    double k_invert = 0.05;     // base coupling magnitude for inversion
    bool normalize_inputs = false; // divide by input count when true
    bool hard_logic = true;     // stronger coupling intent
};

inline double coupling_gain(std::size_t inputs, const CouplingConfig& cfg) {
    if (!cfg.normalize_inputs || inputs == 0) return 1.0;
    return 1.0 / static_cast<double>(inputs);
}

// Result of a single run.
struct PhysResult {
    PhysState final_state;
    uint64_t used_seed = 0;
    int steps = 0;
};

// Result of a coupled (N-parametron) run.
struct CoupledResult {
    std::vector<PhysState> final_states;
    uint64_t used_seed = 0;
    int steps = 0;
};

// Evaluate a single parametron for the given params.
// If emit is true, outputs time series to stdout (CSV/text depending on flags).
PhysResult run_single(const PhysParams& p, bool emit = false);

// Evaluate N coupled parametrons with coupling matrix k (u_i += Σ k_ij * x_j).
CoupledResult run_coupled(const std::vector<PhysParams>& ps, const std::vector<std::vector<double>>& k, bool emit = false);

// Convert physical state to logical bit (0/1).
// Phase-based: 0..π maps to 1, π..2π maps to 0 (normalized to [-π, π]).
inline int bit_from_state(const PhysState& s) {
    constexpr double pi = 3.14159265358979323846;
    double ph = s.phase;
    while (ph <= -pi) ph += 2.0 * pi;
    while (ph > pi) ph -= 2.0 * pi;
    return (ph >= -0.5 * pi && ph <= 0.5 * pi) ? 1 : 0;
}

// Check validity by amplitude threshold.
inline bool valid_from_state(const PhysState& s, double amp_min) {
    return (s.amp >= amp_min);
}

} // namespace parametron_phys
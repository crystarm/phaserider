# Phaserider

A set of research utilities in C++ for modeling and verifying parametron circuits (majority/inverter gates with delays), plus numerical simulation of parametric oscillators. The project is built with CMake and contains several standalone CLI tools.

## Project Contents

The build produces the following binaries (see `CMakeLists.txt`):

- `parametron_ode` — numerical integration of ODEs for two coupled parametric oscillators.
- `parametron_logic` — reference 8‑bit logic adder implementation (majority/inverter) + self‑test.
- `parametron_tick` — tick‑based model of an 8‑bit parametron adder (with delays and tracing).
- `phazeride_adc` — universal N‑bit adder with flags (`C`, `V`, `Z`, `N`) using either the logic model or the tick simulator.
- `circuit_tick` — library (internal core) used by `parametron_tick` and `phazeride_adc`.

Also included:
- `cli.sh` — interactive REPL wrapper around `phazeride_adc` and `parametron_tick`.
- `test.sh` — script comparing `phazeride_adc` against a C reference (exhaustive for 8‑bit, randomized for other widths).
- GitHub Actions: build (`.github/workflows/build.yml`) and CI with tests (`.github/workflows/test.yml`).

## Quick Start

### Dependencies
- CMake ≥ 3.16
- C++20 compiler
- `bash` for `cli.sh` and `test.sh`

### Build

~~~sh
cmake -S . -B build
cmake --build build -j
~~~

Binaries will appear in `build/`.

### Run the CLI REPL

~~~sh
./cli.sh
~~~

The REPL supports adder runs, tracing, and simulation settings (see below).

### Run Tests

~~~sh
./test.sh
~~~

The script:
- builds the project,
- compiles a C reference,
- checks all combinations for width 8,
- runs random tests across selected widths.

Test parameters:

~~~sh
./test.sh --rand 5000 --seed 123456789 --widths "1 2 4 8 16 32 64"
~~~

## Tools

### 1) `phazeride_adc`

An N‑bit adder with flags. It can:
- compute `a + b + cin`,
- compute flags `C`, `V`, `Z`, `N`,
- run in two models: `logic` (combinational) and `tick` (delay simulation),
- print in `dec|hex|bin|all`.

Example:

~~~sh
./build/phazeride_adc --a 250 --b 13 --cin 1 --width 8 --model tick --trace 1
~~~

Key parameters:
- `--a`, `--b` — inputs (decimal, hex `0x..`, bin `0b..`, negatives supported).
- `--cin` — input carry (0/1).
- `--width` — width (1..64).
- `--model` — `logic` or `tick`.
- `--signed` — affects decimal printing only.
- `--fmt` — `dec|hex|bin|all`.
- `--maj_delay`, `--not_delay`, `--ticks`, `--p_flip`, `--seed`, `--trace` — tick simulator parameters.
- `--strict 1` — exit with error if result differs from reference.

Output includes:
- input values,
- sum,
- flags `C V Z N`,
- reference `ref_sum` and `ref_cout`.

### 2) `parametron_tick`

Tick‑based simulation of an 8‑bit parametron adder with delays and tracing. Includes a full test mode.

Examples:

~~~sh
./build/parametron_tick --a 13 --b 250 --cin 1
./build/parametron_tick --a 13 --b 250 --cin 1 --maj_delay 2 --not_delay 1 --ticks 40
./build/parametron_tick --test --maj_delay 1 --not_delay 1
~~~

Parameters:
- `--ticks` — number of ticks (0 = critical).
- `--maj_delay`, `--not_delay` — element delays.
- `--p_flip` — flip‑error probability.
- `--seed` — RNG seed.
- `--trace` — enable trace output.

### 3) `parametron_logic`

Logic (combinational) 8‑bit adder implementation. Useful as a reference.

Examples:

~~~sh
./build/parametron_logic --a 13 --b 250 --cin 1
./build/parametron_logic --test
~~~

### 4) `parametron_ode`

Numerical ODE solver for two coupled parametric oscillators (RK4). Modes:

- `run` — standard run with phase/amplitude tracking.
- `wave` — outputs raw sinusoidal signals.
- `sweep` — sweeps over `h` range to evaluate bistability.

Examples:

~~~sh
./build/parametron_ode --mode run --k_1to2 0.08 --k_2to1 0 --handoff_t 0.05 --k21_off 0 --k12_off 0 --kick_amp2 0 --kick_dur2 0
./build/parametron_ode --mode wave --t 0.001 --dt 1e-6 --out_every 1 --csv 1 > wave.csv
./build/parametron_ode --mode sweep --k 0.05 --h_from 0.2 --h_to 1.1 --h_step 0.02
~~~

Key parameters:
- `--f0`, `--dt`, `--t`, `--zeta`, `--h`, `--alpha_ratio`, `--det_tau`, `--noise`
- `--k12`, `--k21` or `--k` (set both couplings)
- `--handoff_t`, `--k12_off`, `--k21_off`
- `--kick_*` (shape and durations of startup “kicks”)
- `--out_every`, `--csv`
- `--h_from`, `--h_to`, `--h_step`, `--amp_min`
- `--seed`

## `cli.sh` — Interactive REPL

`cli.sh` is a convenient shell for `phazeride_adc` and `parametron_tick`.

Supports:
- `add`, `trace` commands,
- width/model/format and simulation parameters,
- automatic builds.

Example session:

~~~text
phaserider repl. type 'help'.
> width 8
> model tick
> add 13 250 1
> trace 13 250 1
> exit
~~~

Supported settings:
- `width 1..64`
- `model tick|logic`
- `fmt dec|hex|bin|all`
- `signed auto|0|1`
- `maj_delay`, `not_delay`, `ticks`, `p_flip`, `seed`
- `compare on|off` (for `tick` model: compare vs logic)
- `sat off|unsigned|signed` (for `adc`)
- `json 0|1` (JSON‑only output for `adc`)

## Architecture and Implementation

### Adder logic

In `logic.cpp`, basic gates are built from majority and inverter:

- `maj3(a,b,c)` → 1 if majority of inputs are 1.
- `inv1(x)` → inversion.
- `and2`, `or2`, `xor2` are composed from `maj3` and `inv1`.

The adder is a ripple‑carry chain of full adders (8‑bit).

### Tick model (`circuit_tick`)

The “ticker” builds a graph of nodes with delays:
- `kind=0` — constant,
- `kind=1` — input,
- `kind=2` — majority gate,
- `kind=3` — inverter.

Each node has a pipeline of length `delay`, which is the discrete delay. Simulation proceeds by ticks:
- compute next values,
- shift pipelines.

`circuit_tick_crit_ticks(...)` computes the minimal number of ticks required for output stabilization.

### `phazeride_adc`

Combines:
- logic adder (`model=logic`),
- tick simulator (`model=tick`),
- flag computation `C V Z N`,
- comparison against arithmetic reference.

## CI

GitHub Actions:

- `build.yml` — build only.
- `test.yml` — build + run `test.sh`.

## Project Structure

~~~text
.
├── CMakeLists.txt
├── cli.sh
├── test.sh
├── src/
│   ├── adc.cpp
│   ├── circuit_tick.cpp
│   ├── circuit_tick.h
│   ├── logic.cpp
│   ├── ode.cpp
│   └── tick.cpp
└── .github/workflows/
    ├── build.yml
    └── test.yml
~~~

## Notes

- `phazeride_adc` and `cli.sh` support negative numbers and multiple formats (dec/hex/bin).
- Adder width is limited to 64 bits.
- For the tick model, `maj_delay`, `not_delay`, and number of ticks are important; if the result doesn’t converge, increase `ticks` or use `p_flip=0`.

## License

MIT
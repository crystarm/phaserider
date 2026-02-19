#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="$root/build"
adc_bin="$build_dir/phazeride_adc"

widths=(1 2 3 4 5 6 7 8 9 10 12 16 24 32 48 64)
rand_cases=5000
seed=123456789
width8_progress_step=32
progress_step=200

usage() {
  echo "Usage: $0 [--rand N] [--seed N] [--widths \"1 2 4 8 16 32 64\"]"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rand)
      rand_cases="$2"
      shift 2
      ;;
    --seed)
      seed="$2"
      shift 2
      ;;
    --widths)
      IFS=' ' read -r -a widths <<< "$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown arg: $1" >&2
      usage
      exit 1
      ;;
  esac
done

# Build using CMake
mkdir -p "$build_dir"
cmake -S "$root" -B "$build_dir" >/dev/null
cmake --build "$build_dir" -j >/dev/null

if [[ ! -x "$adc_bin" ]]; then
  echo "Missing phazeride_adc after build" >&2
  exit 2
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT
c_ref="$tmpdir/c_ref"

# Compile C reference (NOT C++)
cat > "$tmpdir/c_ref.c" <<'CREF'
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>

static uint64_t mask_for_width(int w) {
    if (w >= 64) return ~0ull;
    return (1ull << w) - 1ull;
}

int main(int argc, char** argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: c_ref width a b cin\n");
        return 2;
    }
    int width = atoi(argv[1]);
    int64_t a_in = (int64_t)strtoll(argv[2], NULL, 10);
    int64_t b_in = (int64_t)strtoll(argv[3], NULL, 10);
    int cin = atoi(argv[4]) & 1;

    uint64_t mask = mask_for_width(width);
    uint64_t a = ((uint64_t)a_in) & mask;
    uint64_t b = ((uint64_t)b_in) & mask;

    __int128 full = ( (__int128)(int64_t)a_in + (__int128)(int64_t)b_in + (cin & 1) );
    uint64_t sum = (uint64_t)full & mask;

    uint64_t sign = (width >= 64) ? (1ull << 63) : (1ull << (width - 1));
    int C = (( (uint64_t)a + (uint64_t)b + (uint64_t)(cin & 1) ) >> width) & 1;
    int V = ((~(a ^ b) & (a ^ sum) & sign) != 0) ? 1 : 0;
    int Z = (sum == 0) ? 1 : 0;
    int N = ((sum & sign) != 0) ? 1 : 0;

    printf("sum=0x%llx C=%d V=%d Z=%d N=%d\n",
           (unsigned long long)sum, C, V, Z, N);
    return 0;
}
CREF

cc -std=c11 -O2 "$tmpdir/c_ref.c" -o "$c_ref"

# Helper to run adc directly and extract sum + flags
run_adc() {
  local width="$1" a="$2" b="$3" cin="$4"

  local out
  out="$("$adc_bin" \
    --a "$a" --b "$b" --cin "$cin" \
    --width "$width" --model tick \
    --fmt hex --signed 1 \
    --maj_delay 1 --not_delay 1 --ticks 0 \
    --p_flip 0.0 --seed 1 --trace 0 \
  )"

  local sum_hex
  local flags
  sum_hex="$(echo "$out" | grep -E '^sum=' | head -n1 | sed -E 's/^sum=0x([0-9a-f]+).*/\1/i')"
  sum_hex="$(echo "$sum_hex" | tr '[:upper:]' '[:lower:]' | sed -E 's/^0+//; s/^$/0/')"
  flags="$(echo "$out" | grep -E '^C=' | head -n1)"

  if [[ -z "$sum_hex" || -z "$flags" ]]; then
    echo "parse_error"
    echo "$out" >&2
    return 3
  fi

  local c v z n
  c="$(echo "$flags" | sed -E 's/.*C=([01]).*/\1/')"
  v="$(echo "$flags" | sed -E 's/.*V=([01]).*/\1/')"
  z="$(echo "$flags" | sed -E 's/.*Z=([01]).*/\1/')"
  n="$(echo "$flags" | sed -E 's/.*N=([01]).*/\1/')"

  echo "sum=0x${sum_hex} C=${c} V=${v} Z=${z} N=${n}"
}

# Compare helper
check_case() {
  local width="$1" a="$2" b="$3" cin="$4"

  local ref
  ref="$("$c_ref" "$width" "$a" "$b" "$cin")"

  local got
  got="$(run_adc "$width" "$a" "$b" "$cin")"

  if [[ "$ref" != "$got" ]]; then
    echo "Mismatch width=$width a=$a b=$b cin=$cin"
    echo "  ref: $ref"
    echo "  got: $got"
    return 1
  fi
  return 0
}

# Deterministic PRNG for random cases
rand_state=$seed
rand64() {
  # xorshift64*
  rand_state=$(( (rand_state ^ (rand_state << 13)) & 0xFFFFFFFFFFFFFFFF ))
  rand_state=$(( (rand_state ^ (rand_state >> 7)) & 0xFFFFFFFFFFFFFFFF ))
  rand_state=$(( (rand_state ^ (rand_state << 17)) & 0xFFFFFFFFFFFFFFFF ))
  echo "$rand_state"
}

echo "Running exhaustive width=8..."
for a in $(seq -128 127); do
  if (( ((a + 128) % width8_progress_step) == 0 )); then
    echo "progress width=8 a=$a"
  fi
  for b in $(seq -128 127); do
    for cin in 0 1; do
      check_case 8 "$a" "$b" "$cin" || exit 1
    done
  done
done
echo "OK width=8 exhaustive"

echo "Running random tests..."
for width in "${widths[@]}"; do
  if [[ "$width" -eq 8 ]]; then
    continue
  fi

  step=$progress_step
  if (( step < 1 )); then step=1; fi
  if (( step > rand_cases )); then step=$rand_cases; fi

  for ((i=0; i<rand_cases; i++)); do
    r1="$(rand64)"
    r2="$(rand64)"
    r3="$(rand64)"

    # map to signed range for width
    # signed min = -2^(w-1), max = 2^(w-1)-1
    if [[ "$width" -ge 63 ]]; then
      # avoid overflow in bash shift
      max=$(( (1<<62) - 1 ))
      min=$(( -1<<62 ))
    else
      max=$(( (1<<(width-1)) - 1 ))
      min=$(( -1<<(width-1) ))
    fi

    a=$(( (r1 % (max - min + 1)) + min ))
    b=$(( (r2 % (max - min + 1)) + min ))
    cin=$(( r3 & 1 ))

    check_case "$width" "$a" "$b" "$cin" || exit 1
    if (( (i + 1) % step == 0 )); then
      echo "progress width=$width $((i + 1))/$rand_cases"
    fi
  done
  echo "OK width=$width random=$rand_cases"
done

echo "ALL TESTS PASSED"
exit 0

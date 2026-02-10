#pragma once

#include <cstdint>

typedef unsigned long long ull;

struct circuit_tick_params
{
    int width;
    int cin;

    int ticks;
    int maj_delay;
    int not_delay;

    double p_flip;
    long long seed;
};

struct circuit_tick_result
{
    ull sum;
    int carry;

    ull used_seed;
    int used_ticks;
    int crit_ticks;
};

typedef void (*circuit_tick_trace_fn)(void* ctx, int t, ull sum, ull carry_pack, int coutv);

int circuit_tick_crit_ticks(int width, int maj_delay, int not_delay);

circuit_tick_result circuit_tick_add(ull a, ull b, const circuit_tick_params& p, circuit_tick_trace_fn trace, void* trace_ctx);

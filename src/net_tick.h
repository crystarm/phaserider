#pragma once

#include <cstdint>

typedef unsigned long long ull;

struct net_tick_params
{
    int width;
    int cin;

    int ticks;
    int maj_delay;
    int not_delay;

    double p_flip;
    long long seed;
};

struct net_tick_result
{
    ull sum;
    int carry;

    ull used_seed;
    int used_ticks;
    int crit_ticks;
};

typedef void (*net_tick_trace_fn)(void* ctx, int t, ull sum, int coutv);

int net_tick_crit_ticks(int width, int maj_delay, int not_delay);

net_tick_result net_tick_add(ull a, ull b, const net_tick_params& p, net_tick_trace_fn trace, void* trace_ctx);

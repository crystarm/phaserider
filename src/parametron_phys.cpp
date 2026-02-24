#include <iostream>
#include <string>
#include <random>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <functional>
#include <vector>

#include "parametron_phys.h"

using namespace std;

typedef long long ll;



struct params
{
    string mode;     // run|wave
    double f0;       // base frequency (Hz)
    double dt;       // time step (s)
    double t_end;    // end time (s)

    double zeta;     // damping
    double h;        // pump strength
    double alpha_ratio; // nonlinear stiffness ratio
    double det_tau;  // detector time constant

    double noise;    // velocity noise level

    std::function<double(double t, double x, double v)> drive_fn;

    int kick_bit;    // initial bit (0/1) encoded as phase
    double kick_amp; // kick amplitude
    double kick_phase; // phase shift for kick (rad)
    double kick_dur; // kick duration

    int out_every;   // output decimation
    int csv;         // csv output

    double amp_min;  // amplitude threshold for valid bit

    int min_periods; // stable periods required
    double amp_eps;  // amplitude change tolerance
    int max_periods; // safety cap

    int seed;        // RNG seed (-1 for time-based)
};









static params default_params()
{
    params p;
    p.mode = "run";
    p.f0 = 50000.0;
    p.dt = 1e-6;
    p.t_end = 0.2;
    p.zeta = 0.010;
    p.h = 0.60;
    p.alpha_ratio = 1.0;
    p.det_tau = 0.010;
    p.noise = 0.0;

    p.kick_bit = 1;
    p.kick_amp = 0.20;
    p.kick_phase = 0.0;
    p.kick_dur = 0.005;

    p.out_every = 10;
    p.csv = 1;

    p.amp_min = 1e-3;

    p.min_periods = 3;
    p.amp_eps = 1e-4;
    p.max_periods = 200;

    p.seed = -1;
    return p;
}

struct deriv_out
{
    double dx;
    double dv;
};

static double kick_drive(double t, double x, double v, const params& p, double w0)
{
    double u = 0.0;
    if (t < p.kick_dur)
    {
        double s = (p.kick_bit != 0) ? +1.0 : -1.0;
        u += s * p.kick_amp * cos(w0 * t + p.kick_phase);
    }
    if (p.drive_fn)
    {
        u += p.drive_fn(t, x, v);
    }
    return u;
}

static double kick_drive_phys(double t, double x, double v, const parametron_phys::PhysParams& p, double w0)
{
    double u = 0.0;
    if (t < p.kick_dur)
    {
        double s = (p.kick_bit != 0) ? +1.0 : -1.0;
        u += s * p.kick_amp * cos(w0 * t + p.kick_phase);
    }
    if (p.drive_fn)
    {
        u += p.drive_fn(t, x, v);
    }
    return u;
}

static deriv_out f_deriv(double t, double x, double v, const params& p, double w0)
{
    double pump = cos(2.0 * w0 * t);
    double alpha = p.alpha_ratio * (w0 * w0);

    double u = kick_drive(t, x, v, p, w0);

    double dx = v;
    double dv = -2.0 * p.zeta * w0 * v
                - (w0 * w0) * (1.0 + p.h * pump) * x
                - alpha * x * x * x
                + u;

    return deriv_out{dx, dv};
}

static void rk4_step(double t, double& x, double& v, const params& p, double w0)
{
    double dt = p.dt;

    auto k1 = f_deriv(t, x, v, p, w0);
    auto k2 = f_deriv(t + 0.5 * dt,
                      x + 0.5 * dt * k1.dx, v + 0.5 * dt * k1.dv,
                      p, w0);
    auto k3 = f_deriv(t + 0.5 * dt,
                      x + 0.5 * dt * k2.dx, v + 0.5 * dt * k2.dv,
                      p, w0);
    auto k4 = f_deriv(t + dt,
                      x + dt * k3.dx, v + dt * k3.dv,
                      p, w0);

    x += (dt / 6.0) * (k1.dx + 2.0 * k2.dx + 2.0 * k3.dx + k4.dx);
    v += (dt / 6.0) * (k1.dv + 2.0 * k2.dv + 2.0 * k3.dv + k4.dv);
}

struct meas
{
    double i;
    double q;
    double amp;
    double phase;
    int bit;
};

struct sim_result
{
    meas m;
    unsigned long long seed;
    double x;
    double v;
};

static sim_result run_sim(const params& p_in, bool emit)
{
    params p = p_in;

    const double pi = acos(-1.0);
    double w0 = 2.0 * pi * p.f0;

    unsigned long long seed = 0;
    if (p.seed >= 0) seed = (unsigned long long)p.seed;
    else seed = (unsigned long long)chrono::high_resolution_clock::now().time_since_epoch().count();

    mt19937_64 rng(seed);
    normal_distribution<double> nd(0.0, 1.0);

    double x = 0.001 * nd(rng);
    double v = 0.001 * nd(rng);

    double beta = p.dt / p.det_tau;
    if (beta > 1.0) beta = 1.0;

    double i = 0.0, q = 0.0;

    ll steps = (ll)llround(p.t_end / p.dt);
    if (steps < 1) steps = 1;

    int steps_per_period = (int)llround(1.0 / (p.f0 * p.dt));
    if (steps_per_period < 1) steps_per_period = 1;

    if (p.max_periods > 0)
    {
        ll max_steps = (ll)p.max_periods * (ll)steps_per_period;
        if (max_steps < steps) steps = max_steps;
    }

    int stable_periods = 0;
    double last_period_amp = 0.0;
    int last_period_bit = 0;
    bool stop = false;

    bool wave = (p.mode == "wave");

    if (emit)
    {
        double spp = 1.0 / (p.f0 * p.dt);
        if (spp < 20.0)
        {
            cerr << "warning: samples_per_period=" << spp << " (<20). dt is likely too large for f0.\n";
        }

        if (p.csv != 0)
        {
            if (wave)
            {
                cout << "t,x,v,pump,drive,i,q,amp,phase,bit\n";
            }
            else
            {
                cout << "t,x,v,pump,i,q,amp,phase,bit\n";
            }
        }
        else
        {
            cout << "seed=" << seed << " f0=" << p.f0 << " dt=" << p.dt << " t=" << p.t_end
                 << " h=" << p.h << " zeta=" << p.zeta
                 << " kick_amp=" << p.kick_amp << " kick_dur=" << p.kick_dur
                 << "\n";
            if (wave)
            {
                cout << "mode=wave (raw) out_every=" << p.out_every << "\n";
            }
        }
    }

    meas last{};
    for (ll k = 0; k <= steps; k++)
    {
        double t = (double)k * p.dt;

        double pump = cos(2.0 * w0 * t);
        double refc = cos(w0 * t);
        double refs = sin(w0 * t);
        double drive = kick_drive(t, x, v, p, w0);

        i = (1.0 - beta) * i + beta * (x * refc);
        q = (1.0 - beta) * q + beta * (x * refs);

        double amp = sqrt(i * i + q * q);
        double phase = atan2(q, i);
        int bit = (phase >= -0.5 * pi && phase <= 0.5 * pi) ? 1 : 0;

        last = meas{i, q, amp, phase, bit};

        if (k > 0 && (k % steps_per_period == 0))
        {
            if (fabs(amp - last_period_amp) <= p.amp_eps && bit == last_period_bit)
            {
                stable_periods++;
            }
            else
            {
                stable_periods = 0;
            }

            last_period_amp = amp;
            last_period_bit = bit;

            if (stable_periods >= p.min_periods) stop = true;
        }

        if (emit && (k % p.out_every == 0 || k == steps))
        {
            if (p.csv != 0)
            {
                if (wave)
                {
                    cout << setprecision(10)
                         << t << "," << x << "," << v << ","
                         << pump << "," << drive << ","
                         << i << "," << q << "," << amp << "," << phase << "," << bit
                         << "\n";
                }
                else
                {
                    cout << setprecision(10)
                         << t << "," << x << "," << v << ","
                         << pump << ","
                         << i << "," << q << "," << amp << "," << phase << "," << bit
                         << "\n";
                }
            }
            else
            {
                if (wave)
                {
                    cout << "t=" << fixed << setprecision(6) << t
                         << " x=" << scientific << setprecision(6) << x
                         << " v=" << scientific << setprecision(6) << v
                         << " pump=" << fixed << setprecision(6) << pump
                         << " drive=" << scientific << setprecision(3) << drive
                         << " bit=" << bit
                         << "\n";
                }
                else
                {
                    cout << "t=" << fixed << setprecision(6) << t
                         << " amp=" << scientific << setprecision(3) << amp
                         << " phase=" << fixed << setprecision(4) << phase
                         << " bit=" << bit
                         << "\n";
                }
            }
        }

        if (stop || k == steps) break;

        rk4_step(t, x, v, p, w0);

        if (p.noise > 0.0)
        {
            double s = p.noise * sqrt(p.dt);
            v += s * nd(rng);
        }
    }

    sim_result r;
    r.m = last;
    r.seed = seed;
    r.x = x;
    r.v = v;
    return r;
}

namespace parametron_phys {



static void compute_deriv_coupled(
    const std::vector<double>& x,
    const std::vector<double>& v,
    const std::vector<PhysParams>& ps,
    const std::vector<std::vector<double>>& k,
    double t,
    double w0,
    std::vector<double>& dx,
    std::vector<double>& dv)
{
    const double pump = cos(2.0 * w0 * t);
    const size_t n = ps.size();

    for (size_t i = 0; i < n; i++)
    {
        double u = kick_drive_phys(t, x[i], v[i], ps[i], w0);
        for (size_t j = 0; j < n; j++)
        {
            u += k[i][j] * x[j];
        }

        const double alpha = ps[i].alpha_ratio * (w0 * w0);
        dx[i] = v[i];
        dv[i] = -2.0 * ps[i].zeta * w0 * v[i]
                - (w0 * w0) * (1.0 + ps[i].h * pump) * x[i]
                - alpha * x[i] * x[i] * x[i]
                + u;
    }
}

CoupledResult run_coupled(const std::vector<PhysParams>& ps, const std::vector<std::vector<double>>& k, bool emit)
{
    CoupledResult out;
    const size_t n = ps.size();
    if (n == 0) return out;
    if (k.size() != n) return out;
    for (size_t i = 0; i < n; i++)
    {
        if (k[i].size() != n) return out;
    }

    const double f0 = ps[0].f0;
    const double dt = ps[0].dt;
    const double w0 = 2.0 * acos(-1.0) * f0;

    for (size_t i = 1; i < n; i++)
    {
        if (ps[i].f0 != f0 || ps[i].dt != dt)
        {
            return out;
        }
    }

    unsigned long long seed = 0;
    if (ps[0].seed >= 0) seed = (unsigned long long)ps[0].seed;
    else seed = (unsigned long long)chrono::high_resolution_clock::now().time_since_epoch().count();

    mt19937_64 rng(seed);
    normal_distribution<double> nd(0.0, 1.0);

    std::vector<double> x(n), v(n), i(n, 0.0), q(n, 0.0);
    std::vector<double> beta(n);
    for (size_t idx = 0; idx < n; idx++)
    {
        x[idx] = 0.001 * nd(rng);
        v[idx] = 0.001 * nd(rng);
        beta[idx] = dt / ps[idx].det_tau;
        if (beta[idx] > 1.0) beta[idx] = 1.0;
    }

    ll steps = (ll)llround(ps[0].t_end / dt);
    if (steps < 1) steps = 1;

    int steps_per_period = (int)llround(1.0 / (f0 * dt));
    if (steps_per_period < 1) steps_per_period = 1;

    int max_periods = ps[0].max_periods;
    for (size_t idx = 1; idx < n; idx++)
    {
        if (ps[idx].max_periods > 0 && (max_periods <= 0 || ps[idx].max_periods < max_periods))
        {
            max_periods = ps[idx].max_periods;
        }
    }
    if (max_periods > 0)
    {
        ll max_steps = (ll)max_periods * (ll)steps_per_period;
        if (max_steps < steps) steps = max_steps;
    }

    std::vector<int> stable_periods(n, 0);
    std::vector<double> last_period_amp(n, 0.0);
    std::vector<int> last_period_bit(n, 0);

    if (emit)
    {
        cout << "t";
        for (size_t idx = 0; idx < n; idx++)
        {
            cout << ",x" << idx << ",v" << idx << ",i" << idx << ",q" << idx
                 << ",amp" << idx << ",phase" << idx << ",bit" << idx;
        }
        cout << "\n";
    }

    std::vector<double> k1x(n), k1v(n), k2x(n), k2v(n), k3x(n), k3v(n), k4x(n), k4v(n);
    std::vector<double> x2(n), v2(n);

    bool stop = false;
    ll kstep = 0;

    for (ll kidx = 0; kidx <= steps; kidx++)
    {
        double t = (double)kidx * dt;

        compute_deriv_coupled(x, v, ps, k, t, w0, k1x, k1v);

        for (size_t idx = 0; idx < n; idx++)
        {
            x2[idx] = x[idx] + 0.5 * dt * k1x[idx];
            v2[idx] = v[idx] + 0.5 * dt * k1v[idx];
        }
        compute_deriv_coupled(x2, v2, ps, k, t + 0.5 * dt, w0, k2x, k2v);

        for (size_t idx = 0; idx < n; idx++)
        {
            x2[idx] = x[idx] + 0.5 * dt * k2x[idx];
            v2[idx] = v[idx] + 0.5 * dt * k2v[idx];
        }
        compute_deriv_coupled(x2, v2, ps, k, t + 0.5 * dt, w0, k3x, k3v);

        for (size_t idx = 0; idx < n; idx++)
        {
            x2[idx] = x[idx] + dt * k3x[idx];
            v2[idx] = v[idx] + dt * k3v[idx];
        }
        compute_deriv_coupled(x2, v2, ps, k, t + dt, w0, k4x, k4v);

        for (size_t idx = 0; idx < n; idx++)
        {
            x[idx] += (dt / 6.0) * (k1x[idx] + 2.0 * k2x[idx] + 2.0 * k3x[idx] + k4x[idx]);
            v[idx] += (dt / 6.0) * (k1v[idx] + 2.0 * k2v[idx] + 2.0 * k3v[idx] + k4v[idx]);
        }

        for (size_t idx = 0; idx < n; idx++)
        {
            if (ps[idx].noise > 0.0)
            {
                double s = ps[idx].noise * sqrt(dt);
                v[idx] += s * nd(rng);
            }

            double refc = cos(w0 * t);
            double refs = sin(w0 * t);
            i[idx] = (1.0 - beta[idx]) * i[idx] + beta[idx] * (x[idx] * refc);
            q[idx] = (1.0 - beta[idx]) * q[idx] + beta[idx] * (x[idx] * refs);
        }

        if (kidx > 0 && (kidx % steps_per_period == 0))
        {
            bool all_stable = true;
            for (size_t idx = 0; idx < n; idx++)
            {
                double amp = sqrt(i[idx] * i[idx] + q[idx] * q[idx]);
                double phase = atan2(q[idx], i[idx]);
                int bit = (phase >= -0.5 * acos(-1.0) && phase <= 0.5 * acos(-1.0)) ? 1 : 0;

                if (fabs(amp - last_period_amp[idx]) <= ps[idx].amp_eps && bit == last_period_bit[idx])
                {
                    stable_periods[idx]++;
                }
                else
                {
                    stable_periods[idx] = 0;
                }

                last_period_amp[idx] = amp;
                last_period_bit[idx] = bit;

                if (stable_periods[idx] < ps[idx].min_periods) all_stable = false;
            }

            if (all_stable) stop = true;
        }

        if (emit && (kidx % ps[0].out_every == 0 || kidx == steps))
        {
            cout << fixed << setprecision(6) << t;
            for (size_t idx = 0; idx < n; idx++)
            {
                double amp = sqrt(i[idx] * i[idx] + q[idx] * q[idx]);
                double phase = atan2(q[idx], i[idx]);
                int bit = (phase >= -0.5 * acos(-1.0) && phase <= 0.5 * acos(-1.0)) ? 1 : 0;
                cout << "," << x[idx] << "," << v[idx] << "," << i[idx] << "," << q[idx]
                     << "," << amp << "," << phase << "," << bit;
            }
            cout << "\n";
        }

        if (stop || kidx == steps)
        {
            kstep = kidx;
            break;
        }
    }

    out.used_seed = seed;
    out.steps = (int)kstep;

    out.final_states.resize(n);
    for (size_t idx = 0; idx < n; idx++)
    {
        double amp = sqrt(i[idx] * i[idx] + q[idx] * q[idx]);
        double phase = atan2(q[idx], i[idx]);
        int bit = (phase >= -0.5 * acos(-1.0) && phase <= 0.5 * acos(-1.0)) ? 1 : 0;

        out.final_states[idx].x = x[idx];
        out.final_states[idx].v = v[idx];
        out.final_states[idx].i = i[idx];
        out.final_states[idx].q = q[idx];
        out.final_states[idx].amp = amp;
        out.final_states[idx].phase = phase;
        out.final_states[idx].bit = bit;
        out.final_states[idx].valid = (amp >= ps[idx].amp_min);
    }

    return out;
}

} // namespace parametron_phys

parametron_phys::PhysResult parametron_phys::run_single(const PhysParams& p, bool emit)
{
    params ip = default_params();
    ip.mode = "run";
    ip.f0 = p.f0;
    ip.dt = p.dt;
    ip.t_end = p.t_end;
    ip.zeta = p.zeta;
    ip.h = p.h;
    ip.alpha_ratio = p.alpha_ratio;
    ip.det_tau = p.det_tau;
    ip.noise = p.noise;
    ip.drive_fn = p.drive_fn;

    ip.kick_bit = p.kick_bit;
    ip.kick_amp = p.kick_amp;
    ip.kick_phase = p.kick_phase;
    ip.kick_dur = p.kick_dur;

    ip.out_every = p.out_every;
    ip.csv = 1;

    ip.amp_min = p.amp_min;
    ip.min_periods = p.min_periods;
    ip.amp_eps = p.amp_eps;
    ip.max_periods = p.max_periods;

    ip.seed = p.seed;

    sim_result r = run_sim(ip, emit);

    PhysResult out;
    out.used_seed = r.seed;
    out.steps = (int)llround(ip.t_end / ip.dt);

    out.final_state.x = r.x;
    out.final_state.v = r.v;
    out.final_state.i = r.m.i;
    out.final_state.q = r.m.q;
    out.final_state.amp = r.m.amp;
    out.final_state.phase = r.m.phase;
    out.final_state.bit = r.m.bit;

    const double amp_min = p.amp_min;
    out.final_state.valid = (out.final_state.amp >= amp_min);

    return out;
}


#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <random>
#include <chrono>
#include <iomanip>
#include <cmath>

using namespace std;

typedef long long ll;

static const int OK = 0;
static const int ERR_USAGE = 1;

struct params
{
    string mode;

    double f0;
    double dt;
    double t_end;

    double zeta;
    double h;

    double alpha_ratio;

    double det_tau;

    double noise;

    int kick_bit1;
    int kick_bit2;
    double kick_amp1;
    double kick_amp2;
    double kick_dur1;
    double kick_dur2;

    double k12;
    double k21;

    double handoff_t;
    double k12_off;
    double k21_off;

    int out_every;
    int csv;

    double h_from;
    double h_to;
    double h_step;

    double amp_min;

    int seed;
};

static int parse_int(const string& v, int& out)
{
    try
    {
        size_t p = 0;
        int x = stoi(v, &p, 10);
        if (p != v.size()) return -1;
        out = x;
        return 0;
    }
    catch (...)
    {
        return -1;
    }
}

static int parse_double(const string& v, double& out)
{
    try
    {
        size_t p = 0;
        double x = stod(v, &p);
        if (p != v.size()) return -1;
        out = x;
        return 0;
    }
    catch (...)
    {
        return -1;
    }
}

static unordered_map<string, string> parse_args(int argc, char** argv)
{
    unordered_map<string, string> m;
    for (int i = 1; i < argc; i++)
    {
        string a = argv[i];
        if (a.rfind("--", 0) == 0)
        {
            string key = a.substr(2);
            string val = "1";
            if (i + 1 < argc)
            {
                string nxt = argv[i + 1];
                if (nxt.rfind("--", 0) != 0)
                {
                    val = nxt;
                    i++;
                }
            }
            m[key] = val;
        }
    }
    return m;
}

static void print_usage()
{
    cout << "usage:\n";
    cout << "  ./parametron_ode --mode run [options]\n";
    cout << "  ./parametron_ode --mode sweep [options]\n\n";

    cout << "common:\n";
    cout << "  --f0 HZ              (default 50000)\n";
    cout << "  --dt SECONDS         (default 1e-6)\n";
    cout << "  --t SECONDS          (default 0.2)\n";
    cout << "  --zeta X             (default 0.010)\n";
    cout << "  --h X                (default 0.60)\n";
    cout << "  --alpha_ratio X      (default 1.0)\n";
    cout << "  --det_tau SECONDS    (default 0.010)\n";
    cout << "  --noise X            (default 0.0)\n";
    cout << "  --seed N             (default: run=-1 time-based, sweep=1+idx)\n\n";

    cout << "coupling:\n";
    cout << "  --k21 X              (default 0.0) 1->2\n";
    cout << "  --k12 X              (default 0.0) 2->1\n";
    cout << "  --k_1to2 X           (alias for k21)\n";
    cout << "  --k_2to1 X           (alias for k12)\n";
    cout << "  --k X                (sets both k12 and k21)\n";
    cout << "  --handoff_t SEC      (default -1 disabled)\n";
    cout << "  --k12_off X          (default 0.0) effective k12 after handoff\n";
    cout << "  --k21_off X          (default 0.0) effective k21 after handoff\n\n";

    cout << "kick:\n";
    cout << "  --kick_bit1 0|1      (default 1)\n";
    cout << "  --kick_bit2 0|1      (default 1)\n";
    cout << "  --kick_amp1 X        (default 0.20)\n";
    cout << "  --kick_amp2 X        (default 0.00)\n";
    cout << "  --kick_dur1 SEC      (default 0.005)\n";
    cout << "  --kick_dur2 SEC      (default 0.000)\n\n";

    cout << "run output:\n";
    cout << "  --out_every N        (default 10)\n";
    cout << "  --csv 0|1            (default 1)\n\n";

    cout << "sweep:\n";
    cout << "  --h_from X           (default 0.0)\n";
    cout << "  --h_to X             (default 1.2)\n";
    cout << "  --h_step X           (default 0.02)\n";
    cout << "  --amp_min X          (default 1e-3)\n\n";

    cout << "examples:\n";
    cout << "  ./parametron_ode --mode run --k_1to2 0.08 --k_2to1 0 --handoff_t 0.05 --k21_off 0 --k12_off 0 --kick_amp2 0 --kick_dur2 0\n";
    cout << "  ./parametron_ode --mode sweep --k 0.05 --h_from 0.2 --h_to 1.1 --h_step 0.02\n";
}

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

    p.kick_bit1 = 1;
    p.kick_bit2 = 1;
    p.kick_amp1 = 0.20;
    p.kick_amp2 = 0.00;
    p.kick_dur1 = 0.005;
    p.kick_dur2 = 0.000;

    p.k12 = 0.0;
    p.k21 = 0.0;

    p.handoff_t = -1.0;
    p.k12_off = 0.0;
    p.k21_off = 0.0;

    p.out_every = 10;
    p.csv = 1;

    p.h_from = 0.0;
    p.h_to = 1.2;
    p.h_step = 0.02;

    p.amp_min = 1e-3;

    p.seed = -1;
    return p;
}

struct deriv_out
{
    double dx1;
    double dv1;
    double dx2;
    double dv2;
};

static void get_eff_k(const params& p, double t, double& k12_eff, double& k21_eff)
{
    k12_eff = p.k12;
    k21_eff = p.k21;
    if (p.handoff_t >= 0.0 && t >= p.handoff_t)
    {
        k12_eff = p.k12_off;
        k21_eff = p.k21_off;
    }
}

static deriv_out f_deriv(double t, double x1, double v1, double x2, double v2, const params& p, double w0)
{
    double pump = cos(2.0 * w0 * t);
    double alpha = p.alpha_ratio * (w0 * w0);

    double u1 = 0.0;
    double u2 = 0.0;

    if (t < p.kick_dur1)
    {
        double s = (p.kick_bit1 != 0) ? +1.0 : -1.0;
        u1 += s * p.kick_amp1 * cos(w0 * t);
    }

    if (t < p.kick_dur2)
    {
        double s = (p.kick_bit2 != 0) ? +1.0 : -1.0;
        u2 += s * p.kick_amp2 * cos(w0 * t);
    }

    double k12_eff = 0.0;
    double k21_eff = 0.0;
    get_eff_k(p, t, k12_eff, k21_eff);

    u1 += k12_eff * x2;
    u2 += k21_eff * x1;

    double dx1 = v1;
    double dv1 = -2.0 * p.zeta * w0 * v1
                 - (w0 * w0) * (1.0 + p.h * pump) * x1
                 - alpha * x1 * x1 * x1
                 + u1;

    double dx2 = v2;
    double dv2 = -2.0 * p.zeta * w0 * v2
                 - (w0 * w0) * (1.0 + p.h * pump) * x2
                 - alpha * x2 * x2 * x2
                 + u2;

    return deriv_out{dx1, dv1, dx2, dv2};
}

static void rk4_step(double t, double& x1, double& v1, double& x2, double& v2, const params& p, double w0)
{
    double dt = p.dt;

    auto k1 = f_deriv(t, x1, v1, x2, v2, p, w0);
    auto k2 = f_deriv(t + 0.5 * dt,
                      x1 + 0.5 * dt * k1.dx1, v1 + 0.5 * dt * k1.dv1,
                      x2 + 0.5 * dt * k1.dx2, v2 + 0.5 * dt * k1.dv2,
                      p, w0);
    auto k3 = f_deriv(t + 0.5 * dt,
                      x1 + 0.5 * dt * k2.dx1, v1 + 0.5 * dt * k2.dv1,
                      x2 + 0.5 * dt * k2.dx2, v2 + 0.5 * dt * k2.dv2,
                      p, w0);
    auto k4 = f_deriv(t + dt,
                      x1 + dt * k3.dx1, v1 + dt * k3.dv1,
                      x2 + dt * k3.dx2, v2 + dt * k3.dv2,
                      p, w0);

    x1 += (dt / 6.0) * (k1.dx1 + 2.0 * k2.dx1 + 2.0 * k3.dx1 + k4.dx1);
    v1 += (dt / 6.0) * (k1.dv1 + 2.0 * k2.dv1 + 2.0 * k3.dv1 + k4.dv1);

    x2 += (dt / 6.0) * (k1.dx2 + 2.0 * k2.dx2 + 2.0 * k3.dx2 + k4.dx2);
    v2 += (dt / 6.0) * (k1.dv2 + 2.0 * k2.dv2 + 2.0 * k3.dv2 + k4.dv2);
}

struct meas
{
    double i1;
    double q1;
    double amp1;
    double phase1;
    int bit1;

    double i2;
    double q2;
    double amp2;
    double phase2;
    int bit2;
};

struct sim_result
{
    meas m;
    unsigned long long seed;
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

    double x1 = 0.001 * nd(rng);
    double v1 = 0.001 * nd(rng);
    double x2 = 0.001 * nd(rng);
    double v2 = 0.001 * nd(rng);

    double beta = p.dt / p.det_tau;
    if (beta > 1.0) beta = 1.0;

    double i1 = 0.0, q1 = 0.0;
    double i2 = 0.0, q2 = 0.0;

    ll steps = (ll)llround(p.t_end / p.dt);
    if (steps < 1) steps = 1;

    if (emit)
    {
        double spp = 1.0 / (p.f0 * p.dt);
        if (spp < 20.0)
        {
            cerr << "warning: samples_per_period=" << spp << " (<20). dt is likely too large for f0.\n";
        }

        if (p.csv != 0)
        {
            cout << "t,x1,v1,x2,v2,pump,refc,refs,k12_eff,k21_eff,i1,q1,amp1,phase1,bit1,i2,q2,amp2,phase2,bit2\n";
        }
        else
        {
            cout << "seed=" << seed << " f0=" << p.f0 << " dt=" << p.dt << " t=" << p.t_end
                 << " h=" << p.h << " zeta=" << p.zeta
                 << " k12=" << p.k12 << " k21=" << p.k21
                 << " handoff_t=" << p.handoff_t
                 << " k12_off=" << p.k12_off << " k21_off=" << p.k21_off
                 << "\n";
        }
    }

    meas last{};
    for (ll k = 0; k <= steps; k++)
    {
        double t = (double)k * p.dt;

        double pump = cos(2.0 * w0 * t);
        double refc = cos(w0 * t);
        double refs = sin(w0 * t);

        i1 = (1.0 - beta) * i1 + beta * (x1 * refc);
        q1 = (1.0 - beta) * q1 + beta * (x1 * refs);

        i2 = (1.0 - beta) * i2 + beta * (x2 * refc);
        q2 = (1.0 - beta) * q2 + beta * (x2 * refs);

        double amp1 = sqrt(i1 * i1 + q1 * q1);
        double phase1 = atan2(q1, i1);
        int bit1 = (i1 >= 0.0) ? 1 : 0;

        double amp2 = sqrt(i2 * i2 + q2 * q2);
        double phase2 = atan2(q2, i2);
        int bit2 = (i2 >= 0.0) ? 1 : 0;

        last = meas{i1, q1, amp1, phase1, bit1, i2, q2, amp2, phase2, bit2};

        if (emit && (k % p.out_every == 0 || k == steps))
        {
            double k12_eff = 0.0;
            double k21_eff = 0.0;
            get_eff_k(p, t, k12_eff, k21_eff);

            if (p.csv != 0)
            {
                cout << setprecision(10)
                     << t << "," << x1 << "," << v1 << "," << x2 << "," << v2 << ","
                     << pump << "," << refc << "," << refs << ","
                     << k12_eff << "," << k21_eff << ","
                     << i1 << "," << q1 << "," << amp1 << "," << phase1 << "," << bit1 << ","
                     << i2 << "," << q2 << "," << amp2 << "," << phase2 << "," << bit2
                     << "\n";
            }
            else
            {
                cout << "t=" << fixed << setprecision(6) << t
                     << " amp1=" << scientific << setprecision(3) << amp1
                     << " bit1=" << bit1
                     << " amp2=" << scientific << setprecision(3) << amp2
                     << " bit2=" << bit2
                     << " k12_eff=" << fixed << setprecision(4) << k12_eff
                     << " k21_eff=" << fixed << setprecision(4) << k21_eff
                     << "\n";
            }
        }

        if (k == steps) break;

        rk4_step(t, x1, v1, x2, v2, p, w0);

        if (p.noise > 0.0)
        {
            double s = p.noise * sqrt(p.dt);
            v1 += s * nd(rng);
            v2 += s * nd(rng);
        }
    }

    sim_result r;
    r.m = last;
    r.seed = seed;
    return r;
}

static int main_run(const params& p)
{
    run_sim(p, true);
    return OK;
}

static int main_sweep(params p)
{
    if (p.h_step <= 0.0) return ERR_USAGE;

    cout << "h,bit1_k0,amp1_k0,bit2_k0,amp2_k0,ok_k0,bit1_k1,amp1_k1,bit2_k1,amp2_k1,ok_k1,bistable1,bistable2\n";

    int base_seed = (p.seed >= 0) ? p.seed : 1;

    int n = (int)floor((p.h_to - p.h_from) / p.h_step + 0.5) + 1;
    if (n < 1) n = 1;

    for (int idx = 0; idx < n; idx++)
    {
        double h = p.h_from + (double)idx * p.h_step;
        if (p.h_step > 0.0 && h > p.h_to + 1e-12) break;

        params p0 = p;
        params p1 = p;
        p0.h = h;
        p1.h = h;

        p0.kick_bit1 = 0;
        p1.kick_bit1 = 1;

        int cur_seed = base_seed + idx;
        p0.seed = cur_seed;
        p1.seed = cur_seed;

        auto r0 = run_sim(p0, false);
        auto r1 = run_sim(p1, false);

        bool ok0_1 = (r0.m.amp1 >= p.amp_min);
        bool ok1_1 = (r1.m.amp1 >= p.amp_min);
        bool ok0_2 = (r0.m.amp2 >= p.amp_min);
        bool ok1_2 = (r1.m.amp2 >= p.amp_min);

        int bistable1 = (ok0_1 && ok1_1 && (r0.m.bit1 != r1.m.bit1)) ? 1 : 0;
        int bistable2 = (ok0_2 && ok1_2 && (r0.m.bit2 != r1.m.bit2)) ? 1 : 0;

        cout << setprecision(10)
             << h << ","
             << r0.m.bit1 << "," << r0.m.amp1 << "," << r0.m.bit2 << "," << r0.m.amp2 << "," << (ok0_1 ? 1 : 0) << ","
             << r1.m.bit1 << "," << r1.m.amp1 << "," << r1.m.bit2 << "," << r1.m.amp2 << "," << (ok1_1 ? 1 : 0) << ","
             << bistable1 << "," << bistable2
             << "\n";
    }

    return OK;
}

int main(int argc, char** argv)
{
    auto args = parse_args(argc, argv);
    params p = default_params();

    if (args.count("mode")) p.mode = args["mode"];

    if (args.count("f0") && parse_double(args["f0"], p.f0) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("dt") && parse_double(args["dt"], p.dt) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("t") && parse_double(args["t"], p.t_end) != 0) { print_usage(); return ERR_USAGE; }

    if (args.count("zeta") && parse_double(args["zeta"], p.zeta) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("h") && parse_double(args["h"], p.h) != 0) { print_usage(); return ERR_USAGE; }

    if (args.count("alpha_ratio") && parse_double(args["alpha_ratio"], p.alpha_ratio) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("det_tau") && parse_double(args["det_tau"], p.det_tau) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("noise") && parse_double(args["noise"], p.noise) != 0) { print_usage(); return ERR_USAGE; }

    if (args.count("kick_bit1") && parse_int(args["kick_bit1"], p.kick_bit1) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("kick_bit2") && parse_int(args["kick_bit2"], p.kick_bit2) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("kick_amp1") && parse_double(args["kick_amp1"], p.kick_amp1) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("kick_amp2") && parse_double(args["kick_amp2"], p.kick_amp2) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("kick_dur1") && parse_double(args["kick_dur1"], p.kick_dur1) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("kick_dur2") && parse_double(args["kick_dur2"], p.kick_dur2) != 0) { print_usage(); return ERR_USAGE; }

    if (args.count("k"))
    {
        double k = 0.0;
        if (parse_double(args["k"], k) != 0) { print_usage(); return ERR_USAGE; }
        p.k12 = k;
        p.k21 = k;
    }

    if (args.count("k_1to2"))
    {
        double k = 0.0;
        if (parse_double(args["k_1to2"], k) != 0) { print_usage(); return ERR_USAGE; }
        p.k21 = k;
    }
    if (args.count("k_2to1"))
    {
        double k = 0.0;
        if (parse_double(args["k_2to1"], k) != 0) { print_usage(); return ERR_USAGE; }
        p.k12 = k;
    }

    if (args.count("k12") && parse_double(args["k12"], p.k12) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("k21") && parse_double(args["k21"], p.k21) != 0) { print_usage(); return ERR_USAGE; }

    if (args.count("handoff_t") && parse_double(args["handoff_t"], p.handoff_t) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("k12_off") && parse_double(args["k12_off"], p.k12_off) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("k21_off") && parse_double(args["k21_off"], p.k21_off) != 0) { print_usage(); return ERR_USAGE; }

    if (args.count("out_every") && parse_int(args["out_every"], p.out_every) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("csv") && parse_int(args["csv"], p.csv) != 0) { print_usage(); return ERR_USAGE; }

    if (args.count("h_from") && parse_double(args["h_from"], p.h_from) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("h_to") && parse_double(args["h_to"], p.h_to) != 0) { print_usage(); return ERR_USAGE; }
    if (args.count("h_step") && parse_double(args["h_step"], p.h_step) != 0) { print_usage(); return ERR_USAGE; }

    if (args.count("amp_min") && parse_double(args["amp_min"], p.amp_min) != 0) { print_usage(); return ERR_USAGE; }

    if (args.count("seed") && parse_int(args["seed"], p.seed) != 0) { print_usage(); return ERR_USAGE; }

    if (p.f0 <= 0.0 || p.dt <= 0.0 || p.t_end <= 0.0 || p.det_tau <= 0.0 || p.out_every <= 0)
    {
        print_usage();
        return ERR_USAGE;
    }

    if (p.mode == "run")
    {
        return main_run(p);
    }
    else if (p.mode == "sweep")
    {
        return main_sweep(p);
    }
    else
    {
        print_usage();
        return ERR_USAGE;
    }
}

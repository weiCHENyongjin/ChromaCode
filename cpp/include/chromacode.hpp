// ChromaCode - C++ reconstruction core (header-only)
// IQ phase-encoded lock-in demodulation for multispectral reconstruction.
//
// A faithful port of the Python `SpectralReconstructor` inverse path:
//   reference -> mix -> zero-phase Butterworth low-pass -> calibrate.
//
// Numerics use Armadillo (vectors / element-wise ops). The Butterworth design
// and the zero-phase forward-backward filter (equivalent to SciPy's
// butter(output='sos') + sosfiltfilt) are implemented here, since no library
// ships a turnkey equivalent.
//
// Dependencies: Armadillo (+ BLAS/LAPACK). Build flags: -larmadillo.
// See ../CODE_STYLE.md for conventions (English identifiers, units in names).

#pragma once

#include <armadillo>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace chromacode {

// ---------------------------------------------------------------------------
// Minimal dependency-free JSON parser (enough for the shared config schema).
// Lets C++ read the SAME .json config Python uses.
// ---------------------------------------------------------------------------
namespace json {

struct Value {
    enum Type { Null, Bool, Num, Str, Arr, Obj } type = Null;
    bool b = false;
    double n = 0.0;
    std::string s;
    std::vector<Value> arr;
    std::map<std::string, Value> obj;

    bool contains(const std::string& k) const {
        return type == Obj && obj.count(k) > 0;
    }
    const Value& at(const std::string& k) const { return obj.at(k); }
    double num() const { return n; }
    double num_or(const std::string& k, double d) const {
        return contains(k) ? obj.at(k).n : d;
    }
};

class Parser {
public:
    explicit Parser(const std::string& text) : s_(text) {}
    Value parse() {
        skip();
        Value v = value();
        return v;
    }

private:
    const std::string& s_;
    std::size_t i_ = 0;

    [[noreturn]] void err(const std::string& m) const {
        throw std::runtime_error("JSON parse error: " + m);
    }
    void skip() {
        while (i_ < s_.size() && std::isspace((unsigned char)s_[i_])) ++i_;
    }
    void expect(const char* lit) {
        std::size_t L = std::strlen(lit);
        if (s_.compare(i_, L, lit) != 0) err(std::string("expected ") + lit);
        i_ += L;
    }
    Value value() {
        skip();
        if (i_ >= s_.size()) err("unexpected end of input");
        char c = s_[i_];
        if (c == '{') return object();
        if (c == '[') return array();
        if (c == '"') { Value v; v.type = Value::Str; v.s = str(); return v; }
        if (c == 't' || c == 'f') return boolean();
        if (c == 'n') { expect("null"); return Value{}; }
        return number();
    }
    Value boolean() {
        Value v; v.type = Value::Bool;
        if (s_[i_] == 't') { expect("true"); v.b = true; }
        else { expect("false"); v.b = false; }
        return v;
    }
    Value number() {
        std::size_t start = i_;
        while (i_ < s_.size() &&
               (std::isdigit((unsigned char)s_[i_]) || std::strchr("+-.eE", s_[i_])))
            ++i_;
        Value v; v.type = Value::Num;
        v.n = std::strtod(s_.substr(start, i_ - start).c_str(), nullptr);
        return v;
    }
    std::string str() {
        ++i_;  // opening quote
        std::string out;
        while (i_ < s_.size() && s_[i_] != '"') {
            char c = s_[i_++];
            if (c == '\\' && i_ < s_.size()) {
                char e = s_[i_++];
                switch (e) {
                    case 'n': out += '\n'; break;
                    case 't': out += '\t'; break;
                    case 'r': out += '\r'; break;
                    case '"': out += '"'; break;
                    case '\\': out += '\\'; break;
                    case '/': out += '/'; break;
                    default: out += e; break;
                }
            } else {
                out += c;
            }
        }
        if (i_ >= s_.size()) err("unterminated string");
        ++i_;  // closing quote
        return out;
    }
    Value array() {
        Value v; v.type = Value::Arr;
        ++i_; skip();
        if (i_ < s_.size() && s_[i_] == ']') { ++i_; return v; }
        while (true) {
            v.arr.push_back(value());
            skip();
            if (i_ >= s_.size()) err("unterminated array");
            if (s_[i_] == ',') { ++i_; continue; }
            if (s_[i_] == ']') { ++i_; break; }
            err("expected ',' or ']'");
        }
        return v;
    }
    Value object() {
        Value v; v.type = Value::Obj;
        ++i_; skip();
        if (i_ < s_.size() && s_[i_] == '}') { ++i_; return v; }
        while (true) {
            skip();
            if (i_ >= s_.size() || s_[i_] != '"') err("expected string key");
            std::string key = str();
            skip();
            if (i_ >= s_.size() || s_[i_] != ':') err("expected ':'");
            ++i_;
            v.obj[key] = value();
            skip();
            if (i_ >= s_.size()) err("unterminated object");
            if (s_[i_] == ',') { ++i_; continue; }
            if (s_[i_] == '}') { ++i_; break; }
            err("expected ',' or '}'");
        }
        return v;
    }
};

inline Value parse(const std::string& text) { return Parser(text).parse(); }
inline Value parse_file(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open config file: " + path);
    std::stringstream ss;
    ss << f.rdbuf();
    return parse(ss.str());
}

}  // namespace json

// ---------------------------------------------------------------------------
// Second-order section (biquad), a0 normalized to 1.
// ---------------------------------------------------------------------------
struct Biquad {
    double b0, b1, b2;  // numerator
    double a1, a2;      // denominator (a0 = 1)
};

// ---------------------------------------------------------------------------
// Butterworth low-pass design -> cascade of biquads (SciPy butter, sos).
//   order : filter order (even; e.g. 4)
//   wn    : cutoff normalized to Nyquist, in (0, 1)
// ---------------------------------------------------------------------------
inline std::vector<Biquad> butterworth_lowpass(int order, double wn) {
    if (order % 2 != 0)
        throw std::invalid_argument("butterworth_lowpass: even order required");
    using cd = std::complex<double>;
    const double fs2 = 4.0;                              // bilinear, fs = 2
    const double warped = 4.0 * std::tan(M_PI * wn / 2.0);

    // Analog prototype poles, scaled (lp2lp); analog gain = warped^order.
    std::vector<cd> pa(order);
    cd prod_4_minus_pa = 1.0;
    for (int k = 0; k < order; ++k) {
        // Butterworth prototype poles (SciPy buttap): m = -N+1, -N+3, ..., N-1.
        const int m = 2 * k - (order - 1);
        cd proto = -std::polar(1.0, M_PI * m / (2.0 * order));
        pa[k] = warped * proto;
        prod_4_minus_pa *= (fs2 - pa[k]);
    }
    double gain = std::pow(warped, order) * std::real(1.0 / prod_4_minus_pa);

    // Bilinear transform: analog poles -> digital poles; zeros all at z = -1.
    std::vector<cd> pd(order);
    for (int k = 0; k < order; ++k)
        pd[k] = (fs2 + pa[k]) / (fs2 - pa[k]);

    // Group conjugate pole pairs into biquads (num zeros (1+z^-1)^2 -> [1,2,1]).
    std::vector<Biquad> sections;
    for (int k = 0; k < order; ++k) {
        if (pd[k].imag() <= 1e-12) continue;            // take one of each conj pair
        double re = pd[k].real();
        double mag2 = std::norm(pd[k]);
        Biquad s{1.0, 2.0, 1.0, -2.0 * re, mag2};
        sections.push_back(s);
    }
    // Put the overall numerator gain into the first section.
    if (!sections.empty()) {
        sections[0].b0 *= gain;
        sections[0].b1 *= gain;
        sections[0].b2 *= gain;
    }
    return sections;
}

// ---------------------------------------------------------------------------
// Zero-phase forward-backward filtering (equivalent to SciPy sosfiltfilt).
// ---------------------------------------------------------------------------
namespace detail {

// Steady-state initial condition of one biquad for a unit-step input.
inline void biquad_zi(const Biquad& s, double& z0, double& z1) {
    const double det = 1.0 + s.a1 + s.a2;
    const double B0 = s.b1 - s.a1 * s.b0;
    const double B1 = s.b2 - s.a2 * s.b0;
    z0 = (B0 + B1) / det;
    z1 = ((1.0 + s.a1) * B1 - s.a2 * B0) / det;
}

// Per-section cascade initial conditions, scaled by cumulative DC gain.
inline std::vector<std::pair<double, double>> cascade_zi(
        const std::vector<Biquad>& sos) {
    std::vector<std::pair<double, double>> zi(sos.size());
    double scale = 1.0;
    for (std::size_t i = 0; i < sos.size(); ++i) {
        double z0, z1;
        biquad_zi(sos[i], z0, z1);
        zi[i] = {scale * z0, scale * z1};
        const double dc = (sos[i].b0 + sos[i].b1 + sos[i].b2)
                        / (1.0 + sos[i].a1 + sos[i].a2);
        scale *= dc;
    }
    return zi;
}

// Filter in place through one biquad (Direct-Form II transposed) with state.
inline void apply_biquad(arma::vec& x, const Biquad& s, double z0, double z1) {
    for (arma::uword i = 0; i < x.n_elem; ++i) {
        const double xi = x[i];
        const double y = s.b0 * xi + z0;
        z0 = s.b1 * xi - s.a1 * y + z1;
        z1 = s.b2 * xi - s.a2 * y;
        x[i] = y;
    }
}

// Odd extension of length `n` at both ends (SciPy padtype='odd').
inline arma::vec odd_ext(const arma::vec& x, int n) {
    const arma::uword N = x.n_elem;
    arma::vec ext(N + 2 * n);
    for (int j = 0; j < n; ++j) ext[j] = 2.0 * x[0] - x[n - j];
    ext.subvec(n, n + N - 1) = x;
    for (int j = 0; j < n; ++j) ext[n + N + j] = 2.0 * x[N - 1] - x[N - 2 - j];
    return ext;
}

}  // namespace detail

inline arma::vec sosfiltfilt(const arma::vec& x, const std::vector<Biquad>& sos) {
    const int nsec = static_cast<int>(sos.size());
    const auto zi = detail::cascade_zi(sos);
    const int padlen = 3 * (2 * nsec + 1);
    if (static_cast<arma::uword>(padlen) >= x.n_elem)
        throw std::invalid_argument("sosfiltfilt: signal too short for padding");

    arma::vec y = detail::odd_ext(x, padlen);

    const double x0f = y[0];
    for (int s = 0; s < nsec; ++s)
        detail::apply_biquad(y, sos[s], zi[s].first * x0f, zi[s].second * x0f);

    y = arma::reverse(y);
    const double x0b = y[0];
    for (int s = 0; s < nsec; ++s)
        detail::apply_biquad(y, sos[s], zi[s].first * x0b, zi[s].second * x0b);
    y = arma::reverse(y);

    return y.subvec(padlen, padlen + x.n_elem - 1);
}

// ---------------------------------------------------------------------------
// Configuration / result types
// ---------------------------------------------------------------------------
struct LEDChannel {
    double wavelength_nm;
    double frequency_hz;
    double phase_rad;        // I-channel = 0, Q-channel = pi/2
    double weight = 1.0;     // per-channel calibration weight
};

struct ReconstructionConfig {
    double sample_rate_hz;
    std::vector<LEDChannel> channels;
    double lpf_cutoff_hz = 3.5;
    int    lpf_order = 4;
    double integration_time_s = -1.0;   // <0 => use integration_delay_s
    double integration_delay_s = -1.0;  // <0 => default (1/fs)/2
};

struct ReconstructionResult {
    std::vector<double> wavelengths_nm;          // channel order
    std::map<double, arma::vec> reflectance;     // wavelength -> series
    arma::vec time_s;
    std::string calibration_mode;
    double integration_delay_s;
};

// Build a config from a parsed JSON object using the SAME schema as Python's
// ReconstructionConfig.from_dict (sensor / lpf / leds). The C++ core uses the
// per-channel `weight`; fields it does not need (e.g. spectral_response) are
// simply ignored, so one config file serves both languages.
inline ReconstructionConfig config_from_json(const json::Value& root) {
    if (root.type != json::Value::Obj || !root.contains("sensor")
            || !root.contains("leds"))
        throw std::runtime_error("config: expected object with 'sensor' and 'leds'");

    ReconstructionConfig cfg;
    const json::Value& sensor = root.at("sensor");
    cfg.sample_rate_hz = sensor.at("fs").num();
    cfg.integration_delay_s =
        sensor.contains("integration_delay") ? sensor.at("integration_delay").num() : -1.0;
    cfg.integration_time_s =
        sensor.contains("integration_time") ? sensor.at("integration_time").num() : -1.0;

    if (root.contains("lpf")) {
        const json::Value& lpf = root.at("lpf");
        cfg.lpf_cutoff_hz = lpf.num_or("cutoff", 3.5);
        cfg.lpf_order = static_cast<int>(lpf.num_or("order", 4));
    }

    for (const json::Value& led : root.at("leds").arr) {
        LEDChannel ch;
        ch.wavelength_nm = led.at("wavelength").num();
        ch.frequency_hz = led.at("freq").num();
        if (led.contains("phase_rad")) ch.phase_rad = led.at("phase_rad").num();
        else if (led.contains("phase_deg")) ch.phase_rad = led.at("phase_deg").num() * M_PI / 180.0;
        else if (led.contains("phase")) ch.phase_rad = led.at("phase").num() * M_PI / 180.0;  // bare = degrees
        else ch.phase_rad = 0.0;
        ch.weight = led.contains("weight") ? led.at("weight").num() : 1.0;
        cfg.channels.push_back(ch);
    }
    return cfg;
}

inline ReconstructionConfig load_config_file(const std::string& path) {
    return config_from_json(json::parse_file(path));
}

// ---------------------------------------------------------------------------
// Reconstructor
// ---------------------------------------------------------------------------
class SpectralReconstructor {
public:
    explicit SpectralReconstructor(ReconstructionConfig config)
        : config_(std::move(config)) {
        const double wn = config_.lpf_cutoff_hz / (config_.sample_rate_hz / 2.0);
        sos_ = butterworth_lowpass(config_.lpf_order, wn);
    }

    // Group delay tau (s) of the sensor's boxcar integration; see theory §6.
    double group_delay_s() const {
        if (config_.integration_delay_s >= 0.0) return config_.integration_delay_s;
        const double exposure = (config_.integration_time_s >= 0.0)
                              ? config_.integration_time_s
                              : 1.0 / config_.sample_rate_hz;
        return 0.5 * exposure;
    }

    // Reconstruct per-wavelength reflectance.
    //   signal          : captured detector series (length N)
    //   white_reference : optional flat-reference capture (empty -> weights mode)
    //   reference_level : known reflectance of the reference target
    //   trim_s          : seconds discarded at each end (IIR edge transient)
    ReconstructionResult reconstruct(const arma::vec& signal,
                                     const arma::vec& white_reference = arma::vec(),
                                     double reference_level = 1.0,
                                     double trim_s = 0.0) const {
        const arma::uword N = signal.n_elem;
        const double fs = config_.sample_rate_hz;
        const arma::vec t = arma::regspace(0, N - 1) / fs;
        const double tau = group_delay_s();

        const bool use_white = white_reference.n_elem == N;
        std::string mode = use_white ? "white_reference" : "weights";

        // Precompute white-reference calibration constants (interior medians).
        std::map<double, double> white_amp;
        if (use_white) {
            const arma::uword lo = static_cast<arma::uword>(0.1 * N);
            const arma::uword hi = static_cast<arma::uword>(0.9 * N) - 1;
            for (const auto& ch : config_.channels) {
                arma::vec ref = reference(ch, t, tau);
                arma::vec s = sosfiltfilt(white_reference % ref, sos_);
                white_amp[ch.wavelength_nm] = arma::median(s.subvec(lo, hi));
            }
        }

        ReconstructionResult result;
        result.time_s = t;
        result.calibration_mode = mode;
        result.integration_delay_s = tau;

        for (const auto& ch : config_.channels) {
            arma::vec ref = reference(ch, t, tau);
            arma::vec lp = sosfiltfilt(signal % ref, sos_);
            arma::vec R;
            if (use_white) {
                R = reference_level * lp / white_amp[ch.wavelength_nm];
            } else {
                R = M_PI * lp / ch.weight;   // R = pi * LPF[s*ref] / w
            }
            result.wavelengths_nm.push_back(ch.wavelength_nm);
            result.reflectance[ch.wavelength_nm] = std::move(R);
        }

        if (trim_s > 0.0) trim(result, trim_s, fs, N);
        return result;
    }

private:
    arma::vec reference(const LEDChannel& ch, const arma::vec& t, double tau) const {
        return arma::sin(2.0 * M_PI * ch.frequency_hz * (t + tau) + ch.phase_rad);
    }

    static void trim(ReconstructionResult& r, double trim_s, double fs,
                     arma::uword N) {
        const arma::uword k = static_cast<arma::uword>(std::round(trim_s * fs));
        if (2 * k >= N) return;
        const arma::uword lo = k, hi = N - 1 - k;
        r.time_s = r.time_s.subvec(lo, hi);
        for (auto& kv : r.reflectance) kv.second = kv.second.subvec(lo, hi);
    }

    ReconstructionConfig config_;
    std::vector<Biquad> sos_;
};

}  // namespace chromacode

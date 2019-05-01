#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

// 0001010
// 0101000

using comp = std::complex<long double>;
constexpr long double pi = 3.14159265359;

unsigned int bit_invert(unsigned int x, int bit_number) {
    x = ((x & 0xffff0000) >> 16) | ((x & 0x0000ffff) << 16);
    x = ((x & 0xff00ff00) >>  8) | ((x & 0x00ff00ff) <<  8);
    x = ((x & 0xf0f0f0f0) >>  4) | ((x & 0x0f0f0f0f) <<  4);
    x = ((x & 0xcccccccc) >>  2) | ((x & 0x33333333) <<  2);
    x = ((x & 0xaaaaaaaa) >>  1) | ((x & 0x55555555) <<  1);

    return x >> (32 - bit_number);
}

std::vector<comp> dft(const std::vector<comp> &f) {
    std::vector<comp> F;

    int N = f.size();
    F.resize(N);

    comp temp;
    for (int t = 0; t < N; t++) {
        F[t] = 0.0;
        for (int x = 0; x < N; x++)
            F[t] += f[x] * std::exp(comp(0.0, -2.0 * pi * x * t / N));
    }

    return F;
}
std::vector<comp> idft(const std::vector<comp> &F) {
    std::vector<comp> f;

    int N = F.size();
    f.resize(N);

    comp temp;
    for (int x = 0; x < N; x++) {
        f[x] = 0.0;
        for (int t = 0; t < N; t++)
            f[x] += F[t] * std::exp(comp(0.0, 2.0 * pi * t * x / N));
        f[x] /= static_cast<long double>(N);
    }

    return f;
}

std::vector<comp> fft(const std::vector<comp> &f) {
    std::vector<comp> F;

    int N = f.size();
    if (N != (N & (-N)))
        return F;

    F.resize(N);
    for (int i = 0; i < N; i++)
        F[i] = f[bit_invert(i, static_cast<int>(std::log2(N)))];

    std::vector<comp> W;
    W.resize(N / 2);
    for (int i = 0; i < W.size(); i++)
        W[i] = std::exp(comp(0.0, 2.0 * pi * i / N));

    comp a, b;
    int M;
    for (int S = 1; S <= N; S *= 2) {
        M = N / (2 * S);

        for (int s = 0; s < S; s++) {
            for (int m = 0; m < M; m++) {
                a = F[2 * S * m + s] + W[M * s] * F[2 * S * m + s + S];
                b = F[2 * S * m + s] - W[M * s] * F[2 * S * m + s + S];

                F[2 * S * m + s    ] = a;
                F[2 * S * m + s + S] = b;
            }
        }
    }

    return F;
}
std::vector<comp> ifft(const std::vector<comp> &f) {
    std::vector<comp> F;

    int N = f.size();
    if (N != (N & (-N)))
        return F;

    F.resize(N);
    for (int i = 0; i < N; i++)
        F[i] = f[bit_invert(i, static_cast<int>(std::log2(N)))];

    std::vector<comp> W;
    W.resize(N / 2);
    for (int i = 0; i < W.size(); i++)
        W[i] = std::exp(comp(0.0, -2.0 * pi * i / N));

    comp a, b;
    int M;
    for (int S = 1; S <= N; S *= 2) {
        M = N / (2 * S);

        for (int s = 0; s < S; s++) {
            for (int m = 0; m < M; m++) {
                a = F[2 * S * m + s] + W[M * s] * F[2 * S * m + s + S];
                b = F[2 * S * m + s] - W[M * s] * F[2 * S * m + s + S];

                F[2 * S * m + s    ] = a;
                F[2 * S * m + s + S] = b;
            }
        }
    }

    for (int i = 0; i < N; i++)
        F[i] = F[i] / static_cast<long double>(N);

    return F;
}

int main() {

    int N = (1 << 15);
    std::vector<comp> f;
    f.resize(N);
    long double freq = 2.0;
    for (int i = 0; i < f.size(); i++) {
        f[i] = 0.0;
        f[i] += 10 * std::cos((1) * 2.0 * pi * static_cast<long double>(i) / static_cast<long double>(N));
        f[i] +=  5 * std::cos((2) * 2.0 * pi * static_cast<long double>(i) / static_cast<long double>(N));
    }

    std::vector<comp> F;
    F = f;

    F = dft(F);
    //F = idft(F);

    //F = fft(F);
    //F = ifft(F);


    for (int i = 0; i < F.size(); i++)
        std::cout << i << " " << std::pow(F[i], 2).real() << std::endl;

    return 0;
}


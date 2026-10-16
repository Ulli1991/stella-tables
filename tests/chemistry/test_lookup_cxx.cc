// C++ port lookup (cxx_lookup_temp_gamma) on the Fortran test points; prints T, gamma.
#include "chem_types.h"
#include <cstdio>
#include <vector>
int main(int argc, char** argv)
{
    std::vector<double> tt(CHEM_ISIZE), gg(CHEM_ISIZE);
    FILE* f = std::fopen(argv[1], "r");
    for (int i = 0; i < CHEM_ISIZE; i++) { double e; if (std::fscanf(f, "%lf %lf %lf", &e, &tt[i], &gg[i]) != 3) return 2; }
    std::fclose(f);
    ChemTables T{}; T.temp_table = tt.data(); T.gamma_table = gg.data();
    f = std::fopen(argv[2], "r"); double a, h, e;
    while (std::fscanf(f, "%lf %lf %lf", &a, &h, &e) == 3) {
        double t, g; cxx_lookup_temp_gamma(a, h, e, &t, &g, T);
        std::printf("%.16e %.16e\n", t, g);
    }
    return 0;
}

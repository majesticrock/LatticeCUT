#pragma once
#include <vector>
#include <cmath>

namespace DOS {
    struct Base {
        virtual ~Base() = default;
        Base(size_t N, double min_energy, double max_energy);

        std::vector<double> const& get_dos();
        inline double get_min_energy() const {
            return _min_energy;
        }
        inline double get_max_energy() const {
            return _max_energy;
        }
    protected:
        static constexpr double _pi = 3.1415926535897932384626;
        virtual void compute() = 0;

        std::vector<double> _dos;
        const double _min_energy;
        const double _max_energy;
        bool _computed{};
    };
}
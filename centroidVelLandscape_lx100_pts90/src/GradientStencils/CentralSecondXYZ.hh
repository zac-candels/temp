#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralSecondXYZ : GradientBase<GradientSecond, Cartesian> {
    template <class TTraits, class TParameter>
    inline double compute(int direction, int k);
};

template <class TTraits, class TParameter>
inline double CentralSecondXYZ::compute(int direction, int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum = 0;

    for (int idx = 1; idx < Stencil::Q; idx++) {
        const double& param1 = TParameter::template get<Lattice>(data.getNeighbor(k, idx));
        const double& param = TParameter::template get<Lattice>(k);

        gradientsum += 2*Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (param1 - param);
    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT * Lattice::DT) * gradientsum;
}

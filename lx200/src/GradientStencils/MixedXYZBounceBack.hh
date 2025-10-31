#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedXYZBounceBack : GradientBase<GradientMixed, Cartesian> {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k);
};

template <class TTraits, class TParameter>
inline double MixedXYZBounceBack::compute(const int direction, const int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    if (this->isBoundary<Lattice>(k)) return 0;

    static DataType& data = DataType::getInstance();

    double gradientsum = 0;
    const static auto& param = TParameter::template get<Lattice>();
    for (int idx = 1; idx < Stencil::Q; idx++) {
        if ((!this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + idx])) &&
            (!this->isBoundary<Lattice>(
                data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx])) &&
            (!this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]))) {
            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 *
                           (-param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx]] +
                            5 * param[data.getNeighbor(k, idx)] - 3 * param[k] -
                            param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]]);

        } else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + idx]))) {
            if (!(this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]])))
                gradientsum +=
                    Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 *
                    (2 * param[k] - 2 * param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]]);

        } else if ((this->isBoundary<Lattice>(
                       data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx]))) {
            if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]))) {
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 *
                               (4 * param[data.getNeighbors()[k * Stencil::Q + idx]] - 4 * param[k]);

            } else
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 *
                               (4 * param[data.getNeighbors()[k * Stencil::Q + idx]] - 3 * param[k] -
                                param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]]);
        } else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]))) {
            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 *
                           (-param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx]] +
                            5 * param[data.getNeighbors()[k * Stencil::Q + idx]] - 4 * param[k]);
        } else
            return 0;
    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;
}

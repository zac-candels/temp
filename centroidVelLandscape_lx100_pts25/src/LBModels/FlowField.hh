#pragma once
#include <utility>

#include "../BoundaryModels/BounceBack.hh"
#include "../Collide.hh"
#include "../Data.hh"
#include "../Forces/ChemicalForce.hh"
#include "../Parameters.hh"
#include "ModelBase.hh"

// FlowField.hh: Contains the details of the LBM model to solve the Navier-Stokes and continuity equation. Each
// Model is given a "TTraits" class that contains stencil, data, force and boundary information

template <class TLattice>
using DefaultTraitFlowField = typename DefaultTrait<TLattice>::template SetBoundary<BounceBack>;

template <class TLattice, class TTraits = DefaultTraitFlowField<TLattice>>
class FlowField : public CollisionBase<TLattice, typename TTraits::Stencil>,
                  public ModelBase<TLattice, TTraits> {  // Inherit from base class to avoid repetition of common
    // calculations

    using Stencil = typename TTraits::Stencil;
    static constexpr int mNDIM = TLattice::NDIM;

   public:
    inline void collide() override;  // Collision step

    inline void initialise() override;  // Initialisation step

    inline void computeMomenta() override;  // Momenta (density, velocity) calculation

    inline double computeEquilibrium(int k, int idx) override;  // Calculate equilibrium in direction idx

    inline void setTau(double val) { mTau = val; mInverseTau = 1.0 / val; }

    template <class, class>
    friend class FlowFieldBinary;

   private:
    double mTau = 1.0;                // TEMPORARY relaxation time
    double mInverseTau = 1.0 / mTau;  // TEMPORARY inverse relaxation time

    std::vector<double>& density = Density<>::get<TLattice>();           // Reference to vector of TDensities
    std::vector<double>& velocity = Velocity<>::get<TLattice, mNDIM>();  // Reference to vector of velocities

    enum { x = 0, y = 1, z = 2 };  // Indices corresponding to x, y, z directions
};

template <class TLattice, class TTraits>
inline void FlowField<TLattice, TTraits>::collide() {  // Collision step

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // loop over k

        if (this->isCollisionNode(k)) {
            double* old_distribution = this->mDistribution.getDistributionOldPointer(k);

            double equilibriums[Stencil::Q];

            for (int idx = 0; idx < Stencil::Q; idx++) {
                equilibriums[idx] = computeEquilibrium(k, idx);
            }

            this->collisionQ(equilibriums, old_distribution, mInverseTau,
                             k);  // CHANGE NEEDED If no forces, don't require them to be passed
        }
    }
}

template <class TLattice, class TTraits>
inline void FlowField<TLattice, TTraits>::initialise() {  // Initialise model
    this->initialiseProcessors();

    this->mData.generateNeighbors();  // Fill array of neighbor values (See Data.hh)
    TTraits::template CollisionModel<Stencil>::template initialise<TLattice>(this->mt_Forces, mTau, mTau);

#pragma omp parallel for schedule(guided)
    for (int k = 0; k < TLattice::N; k++) {  // loop over k

        Density<>::initialise<TLattice>(1.0, k);  // Set density to 1 initially (This will change)
        Velocity<>::initialise<TLattice, mNDIM>(0.0, k, x);
        if constexpr (mNDIM >= 2) Velocity<>::initialise<TLattice, mNDIM>(0.0, k, y);
        if constexpr (mNDIM == 3) Velocity<>::initialise<TLattice, mNDIM>(0.0, k, z);
    }

    ModelBase<TLattice, TTraits>::mData.communicate(Density<>::getInstance<TLattice>());
    ModelBase<TLattice, TTraits>::mData.communicate(Velocity<>::getInstance<TLattice, TTraits::Lattice::NDIM>());

#pragma omp parallel for schedule(guided)
    for (int k = 0; k < TLattice::N; k++) {  // loop over k

        double* distribution = this->mDistribution.getDistributionPointer(k);
        double* old_distribution = this->mDistribution.getDistributionOldPointer(k);

        for (int idx = 0; idx < Stencil::Q; idx++) {
            double equilibrium = computeEquilibrium(k, idx);
            distribution[idx] = equilibrium;  // Set distributions to equillibrium
            old_distribution[idx] = equilibrium;
        }
    }
}

template <class TLattice, class TTraits>
inline void FlowField<TLattice, TTraits>::computeMomenta() {  // Calculate Density<> and Velocity

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // Loop over k

        if (this->isCollisionNode(k)) {
            double* distribution = this->mDistribution.getDistributionPointer(k);
            velocity[k * Stencil::D + x] =
                this->computeVelocity(distribution, this->mt_Forces, density[k], x, k);  // Calculate velocities

            if constexpr (mNDIM >= 2)
                velocity[k * Stencil::D + y] = this->computeVelocity(distribution, this->mt_Forces, density[k], y, k);
            if constexpr (mNDIM == 3)
                velocity[k * Stencil::D + z] = this->computeVelocity(distribution, this->mt_Forces, density[k], z, k);
            density[k] = this->computeDensity(distribution, k);  // Calculate density
        }
    }
}

template <class TLattice, class TTraits>
inline double FlowField<TLattice, TTraits>::computeEquilibrium(int k, int idx) {
    double gamma = CollisionBase<TLattice, Stencil>::computeGamma(&velocity[k * mNDIM], idx);
    return density[k] * gamma;
}

template <class TMethod>
class PressureForce : public ChemicalForceBinaryMu<TMethod> {
   public:
    template <class TTraits>
    inline double computeXYZ(const int xyz, const int k) {
        return TTraits::Stencil::Cs2 *(TTraits::Stencil::Cs2 *
                   GradientDensity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) +
               ChemicalForceBinaryMu<TMethod>::template computeXYZ<TTraits>(xyz, k));
    }
    template <class TTraits>
    inline double computeQ(const int idx, const int k) {
        double sum = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            sum += (TTraits::Stencil::Cs2 *
                        GradientDensity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) +
                    ChemicalForceBinaryMu<TMethod>::template computeXYZ<TTraits>(xyz, k)) *
                   TTraits::Stencil::Ci_xyz(xyz)[idx];
        }
        return TTraits::Stencil::Cs2 *sum;
    }
    template <class TTraits>
    inline double computeDensitySource(int k) {  // SHOULD BE CENTRAL GRADIENTS
        double source = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++)
            source += TTraits::Lattice::DT * 0.5 * TTraits::Stencil::Cs2 *
                      GradientDensity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) *
                      Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz);
        return source;
    }
    template <class TTraitsF>
    inline double computeVelocitySource(const int xyz, const int k) {  // Need to correct velocity

        //return +TTraitsF::Stencil::Cs2*ChemicalPotential<>::get<typename
        //TTraitsF::Lattice>(k)*GradientOrderParameter<>::get<typename
        //TTraitsF::Lattice,TTraitsF::Lattice::NDIM>(k,xyz) * TTraitsF::Lattice::DT / (2.0);
        return -TTraitsF::Stencil::Cs2 * OrderParameter<>::get<typename TTraitsF::Lattice>(k) *
               Gradient<ChemicalPotential<>>::get<typename TTraitsF::Lattice, TTraitsF::Lattice::NDIM>(k, xyz) *
               TTraitsF::Lattice::DT / (2.0);
    }

   private:
};

template <class TMethod>
class PressureForceGamma0 : public ChemicalForceBinaryMu<TMethod> {
   public:
    template <class TTraits>
    inline double computeXYZ(const int xyz, const int k) {
        return -TTraits::Stencil::Cs2 * TTraits::Stencil::Cs2 *
                   GradientDensity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz);
    }
    template <class TTraits>
    inline double computeQ(const int idx, const int k) {
        double sum = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            sum += -(TTraits::Stencil::Cs2 *
                        GradientDensity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz)) *
                   TTraits::Stencil::Ci_xyz(xyz)[idx];
        }
        return TTraits::Stencil::Cs2 * sum;
    }

   private:
};

template <class TLattice, int TNumberOfComponents = 2>
using DefaultTraitFlowFieldPressure =
    typename DefaultTrait<TLattice, TNumberOfComponents>::template SetBoundary<BounceBack>::template SetProcessor<
        Gradients<Density<>, CentralXYZBounceBack>>::template AddForce<PressureForce<He>>;

template <class TLattice, class TTraits = DefaultTraitFlowFieldPressure<TLattice>>
class FlowFieldPressure
    : public CollisionBase<TLattice, typename TTraits::Stencil>,
      public ModelBase<TLattice, TTraits> {  // Inherit from base class to avoid repetition of common
    // calculations

    using Stencil = typename TTraits::Stencil;
    static constexpr int mNDIM = TLattice::NDIM;

   public:

    std::vector<double>& pressure;                  // Reference to vector of TDensities
    std::vector<double>& density;                    // Reference to vector of TDensities
    std::vector<double>& itau;
    std::vector<double>& velocity;  // Reference to vector of velocities

    FlowFieldPressure() : CollisionBase<TLattice,typename TTraits::Stencil>(), ModelBase<TLattice,TTraits>(),
                        pressure(Pressure<>::template get<TLattice>()), density(Density<>::get<TLattice>()),
                        itau(InverseTau<>::get<TLattice>()), velocity(Velocity<>::get<TLattice,TLattice::NDIM>()) {}

    FlowFieldPressure(FlowFieldPressure<TLattice, TTraits>& other) : CollisionBase<TLattice,typename TTraits::Stencil>(), ModelBase<TLattice,TTraits>(other),
                        pressure(Pressure<>::template get<TLattice>()), density(Density<>::get<TLattice>()),
                        itau(InverseTau<>::get<TLattice>()), velocity(Velocity<>::get<TLattice,TLattice::NDIM>()) {}

    FlowFieldPressure<TLattice, TTraits>& operator=(const FlowFieldPressure<TLattice, TTraits>& other) {
        density = Density<>::get<TLattice>();
        pressure = Pressure<>::template get<TLattice>();
        velocity = Velocity<>::get<TLattice,TLattice::NDIM>();
        itau = InverseTau<>::get<TLattice>();
        ModelBase<TLattice,TTraits>::operator=(other);
        return *this;
    }

    inline void setTau1(double val) { mTauMin = val; }
    inline void setTau2(double val) { mTauMax = val; }

    inline void collide() override;  // Collision step

    virtual inline void initialise() override;  // Initialisation step

    virtual inline void computeMomenta() override;  // Momenta (density, velocity) calculation

    virtual inline double computeEquilibrium(int k, int idx) override;  // Calculate equilibrium in direction idx

    enum { x = 0, y = 1, z = 2 };  // Indices corresponding to x, y, z directions

    template <class, int, class>
    friend class FlowFieldPressureNComp;

    template <class, int, class>
    friend class FlowFieldPressureNCompZheng;

    template <class, class>
    friend class FlowFieldPressureWellBalanced;

    template <class, class>
    friend class FlowFieldPressureWellBalanced2;

    template <class, class>
    friend class FlowFieldPressureWellBalanced3Evap;

    template <class, class>
    friend class FlowFieldPressureWellBalanced3;

   private:
    double mTauMin = 1;
    double mTauMax = 1;
};

template <class TLattice, class TTraits>
inline void FlowFieldPressure<TLattice, TTraits>::collide() {  // Collision step

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // loop over k

        if (this->isCollisionNode(k)) {
            double* old_distribution = this->mDistribution.getDistributionOldPointer(k);

            double equilibriums[Stencil::Q];

            for (int idx = 0; idx < Stencil::Q; idx++) {
                equilibriums[idx] = computeEquilibrium(k, idx);
            }

            this->collisionQ(equilibriums, old_distribution, InverseTau<>::get<TLattice>(k),
                             k);  // CHANGE NEEDED If no forces, don't require them to be passed
        }
    }

    ModelBase<TLattice, TTraits>::mData.communicateDistribution();
}

template <class TLattice, class TTraits>
inline void FlowFieldPressure<TLattice, TTraits>::initialise() {  // Initialise model
    this->initialiseProcessors();

    ModelBase<TLattice, TTraits>::mData.generateNeighbors();  // Fill array of neighbor values (See Data.hh)
    TTraits::template CollisionModel<Stencil>::template initialise<TLattice>(this->mt_Forces, mTauMin, mTauMax);

#pragma omp parallel for schedule(guided)
    for (int k = 0; k < TLattice::N; k++) {  // loop over k

        double* distribution = ModelBase<TLattice, TTraits>::mDistribution.getDistributionPointer(k);
        double* old_distribution = ModelBase<TLattice, TTraits>::mDistribution.getDistributionOldPointer(k);
        InverseTau<>::initialise<TLattice>(1.0, k);
        Pressure<>::initialise<TLattice>(1.0/3.0, k);  // Set density to 1 initially (This will change)
        Density<>::initialise<TLattice>(1, k);
        Velocity<>::initialise<TLattice, TLattice::NDIM>(0.0, k, x);
        Velocity<>::initialise<TLattice, TLattice::NDIM>(0.0, k, y);
        if constexpr (TLattice::NDIM == 3) Velocity<>::initialise<TLattice, TLattice::NDIM>(0.0, k, z);

        for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {
            double equilibrium = computeEquilibrium(k, idx);
            distribution[idx] = equilibrium;  // Set distributions to equillibrium
            old_distribution[idx] = equilibrium;
        }
    }

    ModelBase<TLattice, TTraits>::mData.communicate(Pressure<>::getInstance<TLattice>());
    ModelBase<TLattice, TTraits>::mData.communicate(Velocity<>::getInstance<TLattice, TLattice::NDIM>());
}

template <class TLattice, class TTraits>
inline void FlowFieldPressure<TLattice, TTraits>::computeMomenta() {  // Calculate Density<> and Velocity

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // Loop over k

        if (this->isCollisionNode(k)) {
            double* distribution = ModelBase<TLattice, TTraits>::mDistribution.getDistributionPointer(k);

            //pressure[k] = this->computeDensity(distribution, k);  // Calculate density
            
            
            
            // + source; 
            velocity[k * TTraits::Stencil::D + x] =
                //1. / (TTraits::Stencil::Cs2) *
                this->computeVelocity(distribution, this->mt_Forces, density[k], x, k);  // Calculate velocities
            velocity[k * TTraits::Stencil::D + y] =
                //1. / (TTraits::Stencil::Cs2) * 
                this->computeVelocity(distribution, this->mt_Forces, density[k], y, k);
            if constexpr (TLattice::NDIM == 3)
                velocity[k * TTraits::Stencil::D + z] =
                    //1. / (TTraits::Stencil::Cs2) *
                    this->computeVelocity(distribution, this->mt_Forces, density[k], z, k);
            double source = 0;
            for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++)
                source += TTraits::Lattice::DT * 0.5 * 
                        GradientDensity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) *
                        Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz);
            double velfactor = CollisionBase<TLattice, Stencil>::computeVelocityFactor(&velocity[k * mNDIM], 0);
            pressure[k] = (TTraits::Stencil::Cs2/(1.0-Stencil::Weights[0]))*(this->computeDensity(distribution, k)-distribution[0] + source+Stencil::Weights[0]*density[k]*velfactor);


        }
    }

    ModelBase<TLattice, TTraits>::mData.communicate(Pressure<>::getInstance<TLattice>());
    ModelBase<TLattice, TTraits>::mData.communicate(Velocity<>::getInstance<TLattice, TLattice::NDIM>());
}

template <class TLattice, class TTraits>
inline double FlowFieldPressure<TLattice, TTraits>::computeEquilibrium(int k, int idx) {
    double velocityFactor = CollisionBase<TLattice, Stencil>::computeVelocityFactor(&velocity[k * mNDIM], idx);
    return -(idx==0)*pressure[k]/Stencil::Cs2+Stencil::Weights[idx] * (pressure[k]/Stencil::Cs2 + density[k] * velocityFactor);
}
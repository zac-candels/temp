#pragma once
#include <memory>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <utility>

#include "../AddOns/AddOns.hh"
#include "../BoundaryModels/BoundaryBase.hh"
#include "../Collide.hh"
#include "../Data.hh"
#include "../Forces/ForceBase.hh"
#include "../Trait.hh"

/**
 * \file ModelBase.hh This file contains the ModelBase class, which contains useful functions for the LB models
 */

class Model {};  // Used to store the model without template parameters

/**
 * \brief This class contains useful functions for the LBM models, such as collision and boundary computation.
 * \tparam TLattice The lattice class the simulation will run with.
 * \tparam TTraits The traits class for the model.
 * \details This class also stores objects of all the processors and boundary calculators.
 */
template <class TLattice, class TTraits = DefaultTrait<TLattice>>
class ModelBase : public Model {
    static_assert(std::is_base_of<StencilBase, typename TTraits::Stencil>(),
                  "ERROR: invalid TStencil specified in TTraits class.");
    static_assert(
        TTraits::Stencil::D == TLattice::NDIM,
        "ERROR: The chosen TStencil must match the number of TLattice dimensions in the TLattice properties.");
    static_assert(CheckBaseTemplate<ForceBase, typename TTraits::Forces>::value,
                  "ERROR: At least one TForce chosen is not a TForce class. The class must inherit from ForceBase.");

   public:
    ModelBase() : latticeInit(), mData(), mDistribution(mData.getDistributionObject()), mt_Processors(), mt_Forces(), mt_Boundaries() {
        // Initialise the TLattice and parallelisation
    }

    ModelBase(ModelBase<TLattice, TTraits>& other)
        : latticeInit(), mData(other.mData), mDistribution(other.mDistribution), mt_Processors(other.mt_Processors), mt_Forces(other.mt_Forces), mt_Boundaries(other.mt_Boundaries) {}

    ModelBase(ModelBase<TLattice, TTraits>&& other)
        : latticeInit(), mData(other.mData), mDistribution(other.mDistribution), mt_Processors(other.mt_Processors), mt_Forces(other.mt_Forces), mt_Boundaries(other.mt_Boundaries) {}

    ModelBase& operator=(const ModelBase<TLattice, TTraits>& other) {
        mt_Processors=other.mt_Processors;
        mt_Forces=other.mt_Forces;
        mt_Boundaries=other.mt_Boundaries;
        return *this;
    }

    /**
     * \brief Function that will call the compute(k) function for every element of a tuple on each lattice point.
     * \tparam TTupleType The type of the tuple.
     * \param tup Object of the tuple.
     * \param k The lattice index.
     */
    template <class TTupleType>
    inline void computeTuple(TTupleType& tup, int k);

    /**
     * \brief Will compute all the processors for the model over the whole lattice and perform necessary communication
     * steps.
     */
    inline virtual void computeProcessors();  // Perform any necessary computations before collision

    /**
     * brief Will  call runProessor(k) for all elements of the tuple.
     * \tparam TTupleType The type of the tuple.
     * \param tup Object of the tuple.
     * \param k The lattice index.
     */
    template <class TTupleType>
    inline void processorTuple(TTupleType& tup, int k);

    /**
     * brief Will  call communicate() for all elements of the tuple.
     * \tparam TTupleType The type of the tuple.
     * \param tup Object of the tuple.
     */
    template <class TTupleType>
    inline void communicateTuple(TTupleType& tup);

    /**
     * brief Will  call communicate() for all elements of the tuple.
     * \tparam TTupleType The type of the tuple.
     * \param tup Object of the tuple.
     */
    template <class TTupleType>
    inline void communicateBoundaries(TTupleType& tup);

    /**
     * \brief Will  call communicateProcessir() for all elements of the tuple. Intended to be used with tuple of
     * boundary models. \tparam TTupleType The type of the tuple. \param tup Object of the tuple.
     */
    template <class TTupleType>
    inline void communicateProcessorBoundaries(TTupleType& tup);

    /**
     * \brief Pure virtual collide function, must be overriden by LB models. Should perform the entire collision step
     * and possibly the streaming if the stream() function is not used.
     */
    inline virtual void collide() = 0;

    /**
     * \brief Virtual stream function for the streaming step in LBM. May be overriden.
     */
    inline virtual void stream();

    /**
     * \brief Virtual function which applies all the boundary conditions. May be overriden.
     */
    inline virtual void boundaries();  // Boundary calculation

    /**
     * \brief Returns true if any distribution values are nan or infinity. Returns false otherwise.
     */
    inline bool isNan();

    /**
     * \brief Iterates through the lattice and calls the compute function for every boundary calculator in the provided
     * tuple. \tparam TBoundaryType Type of the tuple containing the boundary models. \param boundary Object of the
     * boundary tuple.
     */
    template <class TBoundaryType>
    inline void runBoundaries(TBoundaryType& boundary);

    /**
     * \brief Iterates through the lattice and calls the compute function for every processor or in the provided tuple.
     * \tparam TProcessorType Type of the tuple containing the processors.
     * \param boundary Object of the processor tuple.
     */
    template <class TProcessorType>
    inline void runProcessors(TProcessorType& processor);

    /**
     * \brief Performes any necessary initialisation for the LB model. Must be Overridden.
     */
    inline virtual void initialise() = 0;

    /**
     * \brief Performs the momenta calculation for the LB model (eg density, velocity). Must be overridden.
     */
    inline virtual void computeMomenta() = 0;

    /**
     * \brief Calculates the equilibrium distribution in a given  velocity direction. Must be overridden.
     * \param k The index of the current node in the lattice.
     * \param idx The velocity direction.
     */
    inline virtual double computeEquilibrium(int k, int idx) = 0;

    /**
     * \brief Virtual function to compute the pressure. May be overridden.
     * \param k The index of the current node in the lattice.
     */
    inline virtual double computePressure(int k);

    /**
     * \brief Virtual function to compute the fractional concentration of each fluid (between 0 and 1).
     * \param k The index of the node on the lattice.
     * \param iFluid The index of the fluid, eg. 0 or 1 for binary.
     */
    inline virtual double computeConcentration(int k, int iFluid);

    /**
     * \brief Returns the total force in the chosen cartesian direction on the current lattice node.
     * \param xyz Chosen cartesian direction as an integer (0=x,1=y,2=z).
     * \param k The index of the current node in the lattice.
     */
    inline double computeForces(int xyz, int k) const;

    inline const std::vector<double>& getDistribution()
        const;  //! Returns a reference to the vector containing the distribution functions.

    inline void initialiseProcessors();  //! Calls the initialise() function of the boundaries, forces, and processors.

    /**
     * \brief Returns a unique pointer to a tuple containing objects of all the different forcing methods used by forces
     * in the TTraits class. \tparam TTupleType The type of the tuple containing the forces. \param TForceTuple Tuple
     * containing objects of the forces. \param k The index of the current node in the lattice.
     */
    template <class TTupleType>
    inline auto getForceCalculator(TTupleType& TForceTuple, int k);

    /**
     * \brief Returns an object of the forcing method used by the given force (e.g. Guo forcing).
     * \tparam TForce Type of the force.
     * \param f Reference to object of the force.
     */
    template <class TForce>
    typename TForce::Method getMethod(TForce& f);

    /**
     * \brief Performs any necessary calculations for a force that must take place immediately before the collision
     * loop. \tparam TForce Type of the force. \tparam TForceTuple Type of the tuple containing the forcing methods for
     * all the forces. \param f Reference to object of the force. \param forcemethods Reference to object of the forcing
     * method tuple. \param k The index of the current node in the lattice.
     */
    template <class TForce, typename TForceTuple>
    void precomputeForces(TForce& f, TForceTuple& forcemethods, int k);

    /**
     * \brief Returns a reference to a processor stored in this model with the given type.
     * \tparam TProcessor Type of the processor we want to find.
     * \tparam tuplenum Index of the tuple where the processor lies.
     * \tparam inst Which instance of the processor should the function return (if there are two of the same type).
     * \details This function first checks if the processor is found in the desired tuple, then it calls the
     *          get_type function to create a tuple containing references to all objects of that type in the tuple.
     *          Finally, it returns the desired instance of the processor from that tuple.
     */
    template <class TProcessor, int tuplenum, int inst = 0>
    inline TProcessor& getProcessor();

    /**
     * \brief Returns a reference to the the first processor stored in this model with the given type.
     * \tparam TProcessor Type of the processor we want to find.
     * \details This function first finds the indices of the tuple containing the processor and the
     *          index of the processor in that tuple. It then uses this to find and return the desired
     *          processor.
     */
    template <class TProcessor>
    inline TProcessor& getProcessor();

    /**
     * \brief Returns a reference to a force stored in this model with the given type.
     * \tparam TForce Type of the force we want to find.
     * \tparam inst Which instance of the processor should the function return (if there are two of the same type).
     * \details This function first checks if the processor is found in the TTraits class, then it calls the
     *          get_type function to create a tuple containing references to all objects of that type in the TTraits
     * class. Finally, it returns the desired instance of the processor.
     */
    template <class TForce, int inst = 0>
    inline TForce& getForce();

    /**
     * \brief Returns a reference to the a boundary calculator stored in this model with the given type.
     * \tparam TBoundary Type of the boundary calculator we want to find.
     * \tparam tuplenum Index of the tuple where the boundary calculator lies.
     * \tparam inst Which instance of the boundary calculator should the function return (if there are two of the same
     * type). \details This function first checks if the boundary calculator is found in the desired tuple, then it
     * calls the get_type function to create a tuple containing references to all objects of that type in the tuple.
     *          Finally, it returns the desired instance of the boundary calculator from that tuple.
     */
    template <class TBoundary, int tuplenum, int inst = 0>
    inline TBoundary& getBoundary();

    /**
     * \brief Returns a reference to the the first boundary calculator stored in this model with the given type.
     * \tparam TBoundary Type of the boundary calculator we want to find.
     * \details This function first finds the indices of the tuple containing the boundary calculator and the
     *          index of the boundary calculator in that tuple. It then uses this to find and return the desired
     *          boundary calculator.
     */
    template <class TBoundary>
    inline TBoundary& getBoundary();

    /**
     * \brief This function sums the forces, grouped by the prefactor of the force method that depends on the relaxation
     * time. \tparam maptype Type of the compile time map of prefactor type to an array containing the force sums for
     * each velocity direction. \tparam prefactortuple Type of the tuple containing a struct for each prefactor type
     * with the type and an array to hold the sum \tparam forcetype Type of the force for which we are adding to the sum
     * stored in the prefactors tuple. \param prefactors object of the tuple containing the type and sums for each
     * prefactor. \param f object of the force which will be computed and added to the sum in each direction. \param idx
     * The velocity index. \param k The index of the current node in the lattice. \details We get the struct from the
     * prefactors tuple with the prefactor type matching the given force. Within this struct, there is an array with a
     * number of entries equal to the number of velocity directions. We add to the value stored in these with the
     * computed value of the force in each velocity direction.
     */
    template <class maptype, class prefactortuple, class forcetype>
    inline void setForceSums(prefactortuple& prefactors, forcetype& f, int idx, int k);

    /**
     * \brief This function performs the collision (and potentially streaming) step in a given velocity direction.
     * \param equilibriums Pointer to the first element of an array containing the calculated equilibrium distributions.
     * \param olddistributions Pointer to the first element of an array contaning the distributions from the previous
     * timestep \param inversetau Inverse viscous relaxationn time for this collision. \param k The index of the current
     * node in the lattice. \details This function first saves the equilibrium distributions if needed. Then it
     * determines whether we have any forces. If we do, it creates a unique pointer to a tuple grouping the forces by
     *          their methods. It then creates a compile time map of the viscosity dependent prefactor
     *          for each force to a struct with an array containing the sum of forces for that prefactor.
     */
    inline void collisionQ(const double* equilibriums, const double* olddistributions, const double& inversetau, int k);

    template <class TForceTypes>
    void updateForces(double& TForce, TForceTypes& forcetypes, int k, int idx);

    inline double computeDensity(const double* distribution, const int k);  // Calculate density

    static inline double computeVelocity(const double* distribution, typename TTraits::Forces& forcetuple,
                                         const double& density, const int xyz, const int k);  // Calculate velocity

    TLattice latticeInit;
    typename TTraits::template DataType<TLattice, typename TTraits::Stencil> mData;
    typename TTraits::template DataType<TLattice, typename TTraits::Stencil>::DistributionData& mDistribution =
        mData.getDistributionObject();
    // Distributions

    enum { x = 0, y = 1, z = 2 };  // Indices corresponding to x, y, z directions

    typename TTraits::Processors mt_Processors;
    typename TTraits::Forces mt_Forces;
    typename TTraits::Boundaries mt_Boundaries;
    Geometry<TLattice> mGeometry;

    inline void setCollideID(int id) { mCollideIDs[0] = id; };

    inline void setCollideID(const std::vector<int>& id) { mCollideIDs = id; };

    inline bool isCollisionNode(int k);

    std::vector<int> mCollideIDs = {0};

    std::vector<double>& distribution = mDistribution.getDistribution();  // Reference to vector of distributions
};

template <class TLattice, class TTraits>
template <class TProcessor, int tuplenum, int inst>
inline TProcessor& ModelBase<TLattice, TTraits>::getProcessor() {
    static_assert(has_type<TProcessor, std::tuple_element_t<tuplenum, typename TTraits::Processors>>::value,
                  "Desired processor is not included in the model.");

    auto processors = get_type<TProcessor>(std::get<tuplenum>(mt_Processors));

    return std::get<inst>(processors);
}

template <class TLattice, class TTraits>
template <class TProcessor>
inline TProcessor& ModelBase<TLattice, TTraits>::getProcessor() {
    using indices = typename tuple_tuple_index_of<typename TTraits::Processors, TProcessor>::idx;

    static constexpr int idx1 = std::tuple_element<0, indices>::type::value;
    static constexpr int idx2 = std::tuple_element<1, indices>::type::value;

    return std::get<idx2>(std::get<idx1>(mt_Processors));
}

template <class TLattice, class TTraits>
template <class TForce, int inst>
inline TForce& ModelBase<TLattice, TTraits>::getForce() {
    static_assert(has_type<TForce, typename TTraits::Forces>::value, "Desired force is not included in the model.");

    auto forces = get_type<TForce>(mt_Forces);

    return std::get<inst>(forces);
}

template <class TLattice, class TTraits>
template <class TBoundary, int tuplenum, int inst>
inline TBoundary& ModelBase<TLattice, TTraits>::getBoundary() {
    static_assert(has_type<TBoundary, std::tuple_element_t<tuplenum, typename TTraits::Boundaries>>::value,
                  "Desired boundary is not included in the model.");

    auto boundaries = get_type<TBoundary>(std::get<tuplenum>(mt_Boundaries));

    return std::get<inst>(boundaries);
}

template <class TLattice, class TTraits>
template <class TBoundary>
inline TBoundary& ModelBase<TLattice, TTraits>::getBoundary() {
    using indices = typename tuple_tuple_index_of<typename TTraits::Boundaries, TBoundary>::idx;

    static constexpr int idx1 = std::tuple_element<0, indices>::type::value;
    static constexpr int idx2 = std::tuple_element<1, indices>::type::value;

    return std::get<idx2>(std::get<idx1>(mt_Boundaries));
}

template <class TLattice, class TTraits>
inline void ModelBase<TLattice, TTraits>::initialiseProcessors() {
    // Boundaries
    std::apply(
        [this](auto&... boundarytuple) {
            (std::apply([this](auto&... boundary) { (boundary.initialise(this), ...); }, boundarytuple), ...);
        },
        mt_Boundaries);

    // Processors
    std::apply(
        [this](auto&... processortuple) {
            (std::apply([this](auto&... processor) { (processor.initialise(this), ...); }, processortuple), ...);
        },
        mt_Processors);

    // Forces
    std::apply([this](auto&... force) { (force.initialise(this), ...); }, mt_Forces);
}

template <class TLattice, class TTraits>
inline double ModelBase<TLattice, TTraits>::computeDensity(const double* distribution,
                                                           const int k) {  // Density<> calculation
    // Density<> is the sum of distributions plus any source/correction terms

    if constexpr (std::tuple_size<typename TTraits::Forces>::value != 0) {
        return (CollisionBase<TLattice, typename TTraits::Stencil>::computeZerothMoment(distribution) +
                std::apply([k](auto&... forces) { return (forces.template computeDensitySource<TTraits>(k) + ...); },
                           mt_Forces)) *
               std::apply(
                   [k](auto&... forces) {
                       return (forces.template computeDensitySourceMultiplicative<TTraits>(k) * ...);
                   },
                   mt_Forces);

    } else
        return CollisionBase<TLattice, typename TTraits::Stencil>::computeZerothMoment(distribution);
}

template <class TLattice, class TTraits>
inline double ModelBase<TLattice, TTraits>::computeVelocity(const double* distribution,
                                                            typename TTraits::Forces& forcetuple, const double& density,
                                                            const int xyz,
                                                            const int k) {  // Velocity calculation in direction xyz
    // Velocity in direction xyz is sum of distribution times the xyz component of the discrete velocity vector
    // in each direction plus any source/correction terms

    if constexpr (std::tuple_size<typename TTraits::Forces>::value != 0) {
        return (1. / (Density<>::get<TLattice>(k))) *
                   CollisionBase<TLattice, typename TTraits::Stencil>::computeFirstMoment(distribution, xyz) +
               (1. / (Density<>::get<TLattice>(k))) *
                   std::apply(
                       [xyz, k](auto&&... forces) mutable {
                           return (forces.template computeVelocitySource<TTraits>(xyz, k) + ...);
                       },
                       forcetuple);

    } else
        return (1. / (Density<>::get<TLattice>(k))) *
               CollisionBase<TLattice, typename TTraits::Stencil>::computeFirstMoment(distribution, xyz);
}

template <class TLattice, class TTraits>
inline double ModelBase<TLattice, TTraits>::computePressure(int k) {
    double density = computeDensity(this->mDistribution.getDistributionPointer(k), k);
    return density * TTraits::Stencil::Cs2;
}

template <class TLattice, class TTraits>
inline double ModelBase<TLattice, TTraits>::computeConcentration(int k, int iFluid) {
    return 1;
}

template <class TLattice, class TTraits>
inline const std::vector<double>& ModelBase<TLattice, TTraits>::getDistribution() const {
    return distribution;  // Return reference to distribution vector
}

template <class TLattice, class TTraits>
template <class TTupleType>
inline auto ModelBase<TLattice, TTraits>::getForceCalculator(TTupleType& TForceTuple, int k) {
    if constexpr (std::tuple_size<TTupleType>::value != 0) {
        auto tempforce = std::apply(
            [k, this](auto&... forces) {  // See Algorithm.hh for explanation of std::apply
                return std::make_unique<decltype(make_tuple_unique(std::make_tuple(getMethod(forces))...))>();

            },
            TForceTuple);

        tempforce = std::apply(
            [this, tempforce = std::move(tempforce),
             k](auto&... forces) mutable {  // See Algorithm.hh for explanation of std::apply
                (this->precomputeForces(forces, *tempforce, k), ...);
                return std::move(tempforce);

            },
            TForceTuple);

        return tempforce;

    } else
        return std::make_unique<std::tuple<>>();
}

template <class TLattice, class TTraits>
template <class TForceTypes>
void ModelBase<TLattice, TTraits>::updateForces(double& TForce, TForceTypes& forcetypes, int k, int idx) {
    if constexpr (std::tuple_size<TForceTypes>::value != 0) {
        TForce = std::apply(
            [idx, k](auto&... forcetype) { return (forcetype.template compute<TTraits>(idx, k) + ...); }, forcetypes);
    }
}

template <class TLattice, class TTraits>
template <class maptype, class prefactortuple, class forcetype>
inline void ModelBase<TLattice, TTraits>::setForceSums(prefactortuple& prefactors, forcetype& f, int idx, int k) {
    std::get<typename maptype::template get<typename forcetype::Prefactor>>(prefactors).val[idx] +=
        f.template compute<TTraits>(idx, k);
}

template <class TLattice, class TTraits>
template <class TForce, typename TForceTuple>
void ModelBase<TLattice, TTraits>::precomputeForces(TForce& f, TForceTuple& forcemethods, int k) {
    std::get<decltype(getMethod(f))>(forcemethods).template precompute<TTraits>(f, k);
}

template <class TLattice, class TTraits>
inline void ModelBase<TLattice, TTraits>::stream() {
    TLattice::ResetParallelTracking();

    if constexpr (TTraits::template DataType<TLattice, typename TTraits::Stencil>::IsStreamingSeperate) {
#pragma omp for schedule(guided)
        for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // loop over k

            for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {
                mDistribution.getPostStreamingDistribution(k, idx) = mDistribution.getPostCollisionDistribution(k)[idx];
            }
        }
    }

    this->mData.communicateDistribution();

    std::apply([this](auto&... boundaryprocessor) { (communicateBoundaries(boundaryprocessor), ...); }, mt_Boundaries);
    std::apply([this](auto&... boundaryprocessor) { (communicateTuple(boundaryprocessor), ...); }, mt_Boundaries);

    TLattice::ResetParallelTracking();
}

template <class TLattice, class TTraits>
inline void ModelBase<TLattice, TTraits>::collisionQ(const double* equilibriums, const double* olddistributions,
                                                     const double& inversetau, int k) {
    if (mDistribution.SaveEquilibrium) {
        mDistribution.saveEquilibriums(equilibriums, k);
    }

    if constexpr (std::tuple_size<typename TTraits::Forces>::value != 0) {
        auto forcemethods = getForceCalculator(mt_Forces, k);

        auto tempMap = std::apply(
            [this](auto&... forces) {  // See Algorithm.hh for explanation of std::apply
                ct_map_types<kv<typename std::remove_reference<decltype(forces)>::type::Method::Prefactor,
                                std::array<double, TTraits::Stencil::Q>>...>
                    tempmap;

                return tempmap;

            },
            mt_Forces);

        using ForcingMap = decltype(tempMap);

        auto tempTuple = std::apply(
            [this](auto&... forces) {  // See Algorithm.hh for explanation of std::apply
                constexpr std::tuple<typename ForcingMap::template get<
                    typename std::remove_reference<decltype(forces)>::type::Method::Prefactor>...>
                    temptup;
                constexpr auto temptup2 = make_tuple_unique(temptup);
                return temptup2;

            },
            mt_Forces);
#pragma omp simd
        for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {
            std::apply(
                [this, idx, k, &tempTuple](auto&... forces) mutable {
                    (this->setForceSums<ForcingMap>(tempTuple, forces, idx, k), ...);
                },
                *forcemethods);
        }

        auto tempforceprefactors = std::apply(
            [](auto&... methods) {  // See Algorithm.hh for explanation of std::apply
                return std::make_unique<decltype(make_tuple_unique(std::make_tuple(getForcePrefactor(methods))...))>();

            },
            *forcemethods);
#pragma omp simd
        for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {  // loop over discrete velocity directions
            // Set distribution at location "mDistribution.streamIndex" equal to the value returned by
            //"computeCollisionQ"

            double collision =
                TTraits::template CollisionModel<typename TTraits::Stencil>::template collide<
                    typename TTraits::Lattice>(olddistributions, equilibriums, inversetau, idx) +
                std::apply(
                    [&inversetau, &tempTuple, idx, this](auto&... prefactors) mutable {
                        return (TTraits::template CollisionModel<typename TTraits::Stencil>::template forcing<
                                    typename TTraits::Lattice, decltype(prefactors)>(
                                    this->mt_Forces,
                                    &(std::get<typename ForcingMap::template get<
                                          typename remove_const_and_reference<decltype(prefactors)>::type>>(tempTuple)
                                          .val[0]),
                                    inversetau, idx) +
                                ...);
                    },
                    *tempforceprefactors);

            mDistribution.getPostCollisionDistribution(k, idx) = collision;
        }
    } else {
        for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {  // loop over discrete velocity directions
            // Set distribution at location "mDistribution.streamIndex" equal to the value returned by
            //"computeCollisionQ"

            double collision = TTraits::template CollisionModel<typename TTraits::Stencil>::template collide<
                typename TTraits::Lattice>(olddistributions, equilibriums, inversetau, idx);

            mDistribution.getPostCollisionDistribution(k, idx) = collision;
        }
    }
}

template <class TLattice, class TTraits>
template <class TTupleType>
inline void ModelBase<TLattice, TTraits>::computeTuple(TTupleType& tup, int k) {
    if constexpr (std::tuple_size<TTupleType>::value != 0) {  // Check if there is at least one element
                                                              // in F

        std::apply(
            [k](auto&... obj) {  // See Algorithm.hh for explanation of std::apply
                (obj.template compute<TTraits>(k), ...);

            },
            tup);
    }
}

template <class TLattice, class TTraits>
template <class TTupleType>
inline void ModelBase<TLattice, TTraits>::communicateTuple(TTupleType& tup) {
    if constexpr (std::tuple_size<TTupleType>::value != 0) {  // Check if there is at least one element
                                                              // in F

        std::apply(
            [](auto&... obj) {  // See Algorithm.hh for explanation of std::apply
                (obj.template communicate<TTraits>(), ...);

            },
            tup);
    }
}

template <class TLattice, class TTraits>
template <class TTupleType>
inline void ModelBase<TLattice, TTraits>::communicateBoundaries(TTupleType& tup) {
    if constexpr (std::tuple_size<TTupleType>::value != 0) {  // Check if there is at least one element
                                                              // in F

        std::apply(
            [this](auto&... obj) {  // See Algorithm.hh for explanation of std::apply
                (obj.template communicate<TTraits>(this->mDistribution), ...);

            },
            tup);
    }
}

template <class TLattice, class TTraits>
template <class TTupleType>
inline void ModelBase<TLattice, TTraits>::communicateProcessorBoundaries(TTupleType& tup) {
    if constexpr (std::tuple_size<TTupleType>::value != 0) {  // Check if there is at least one element
                                                              // in F

        std::apply(
            [](auto&... obj) {  // See Algorithm.hh for explanation of std::apply
                (obj.template communicateProcessor<TTraits>(), ...);

            },
            tup);
    }
}

template <class TLattice, class TTraits>
template <class TTupleType>
inline void ModelBase<TLattice, TTraits>::processorTuple(TTupleType& tup, int k) {
    if constexpr (std::tuple_size<TTupleType>::value != 0) {  // Check if there is at least one element
                                                              // in F

        std::apply(
            [k](auto&... obj) {  // See Algorithm.hh for explanation of std::apply
                (obj.template runProcessor<TTraits>(k), ...);

            },
            tup);
    }
}

template <class TLattice, class TTraits>
inline void ModelBase<TLattice, TTraits>::computeProcessors() {
    TLattice::ResetParallelTracking();

    std::apply([this](auto&... processor) { (runProcessors(processor), ...); }, mt_Processors);

    std::apply([this](auto&... boundaryprocessor) { (communicateProcessorBoundaries(boundaryprocessor), ...); },
               mt_Boundaries);

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // loop over k

        std::apply([this, k](auto&... boundaryprocessor) { (processorTuple(boundaryprocessor, k), ...); },
                   mt_Boundaries);

        processorTuple(mt_Forces, k);
    }

    communicateTuple(mt_Forces);

    TLattice::ResetParallelTracking();

#pragma omp master
    {
        mDistribution.getDistribution().swap(mDistribution.getDistributionOld());  // swap old and new distributions
                                                                                   // before collision
    }
#pragma omp barrier
}

template <class TLattice, class TTraits>
inline void ModelBase<TLattice, TTraits>::boundaries() {
    TLattice::ResetParallelTracking();

    std::apply([this](auto&... boundaryprocessor) { (runBoundaries(boundaryprocessor), ...); }, mt_Boundaries);

    TLattice::ResetParallelTracking();
}

template <class TLattice, class TTraits>
template <class TBoundaryType>
inline void ModelBase<TLattice, TTraits>::runBoundaries(TBoundaryType& boundaryprocessor) {
    if constexpr (std::tuple_size<typename TTraits::Boundaries>::value != 0) {
#pragma omp for schedule(guided)
        for (int k = 0; k < TLattice::N; k++) {  // loop over k
            // Check if there are any boundary
            // models
            std::apply(
                [this, k](auto&... boundaries) { (boundaries.template compute<TTraits>(this->mDistribution, k), ...); },
                boundaryprocessor);
        }
    }
}

template <class TLattice, class TTraits>
template <class TProcessorType>
inline void ModelBase<TLattice, TTraits>::runProcessors(TProcessorType& processortuple) {
    if constexpr (std::tuple_size<typename TTraits::Processors>::value != 0) {
#pragma omp for schedule(guided)
        for (int k = 0; k < TLattice::N; k++) {  // loop over k

            std::apply([this, k](auto&... processor) { (processor.template compute<TTraits>(k), ...); },
                       processortuple);
        }
        communicateTuple(processortuple);
    }
}

template <class TLattice, class TTraits>
inline bool ModelBase<TLattice, TTraits>::isNan() {
    bool isnan = false;
#pragma omp for schedule(guided)
    for (auto i : distribution) {
        if (std::isnan(i) || std::isinf(i)) isnan = true;
    }
    return isnan;
}

template <class TLattice, class TTraits>
inline bool ModelBase<TLattice, TTraits>::isCollisionNode(int k) {
    for (int i : mCollideIDs) {
        if (Geometry<TLattice>::getBoundaryType(k) == i) return true;
    }
    return false;
}

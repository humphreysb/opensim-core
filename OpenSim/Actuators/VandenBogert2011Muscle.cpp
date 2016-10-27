/* -------------------------------------------------------------------------- *
 *                      OpenSim:  VandenBogert.cpp                            *
 * -------------------------------------------------------------------------- */


//=============================================================================
// INCLUDES
//=============================================================================
#include "VandenBogert2011Muscle.h"
#include <OpenSim/Simulation/Model/Model.h>


//=============================================================================
// STATICS
//=============================================================================
using namespace std;
using namespace OpenSim;



//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
//_____________________________________________________________________________
/*
 * Default constructor
 */
VandenBogert2011Muscle::VandenBogert2011Muscle()
{
    constructProperties();
    finalizeFromProperties();
}
//_____________________________________________________________________________
/*
 * Constructor.
 */
VandenBogert2011Muscle::VandenBogert2011Muscle
        (const std::string &name,  double maxIsometricForce,
                 double optimalFiberLength,double tendonSlackLength,
                 double pennationAngle)
{ ;
    constructProperties();
    setName(name);
    setMaxIsometricForce(maxIsometricForce);
    setOptimalFiberLength(optimalFiberLength);
    setTendonSlackLength(tendonSlackLength);
    setPennationAngleAtOptimalFiberLength(pennationAngle);
    finalizeFromProperties();
}



//_____________________________________________________________________________
/*
 * Construct and initialize properties.
 * All properties are added to the property set. Once added, they can be
 * read in and written to files.
 */
void VandenBogert2011Muscle::constructProperties()
{

    constructProperty_tendon_strain_at_max_iso_force(0.033);
    constructProperty_active_force_length_curve_width(0.63);
    constructProperty_force_velocity_hill_constant(0.25);
    constructProperty_force_velocity_max_lengthening_force_norm(1.5);
    constructProperty_fiber_damping(0.01);
    constructProperty_fiber_slack_length_norm(1.0);
    constructProperty_activation_time_constant(0.01);
    constructProperty_deactivation_time_constant(0.04);
    constructProperty_pennation_at_optimal_fiber_length(0);
    constructProperty_default_activation(0.1);
    constructProperty_default_fiber_length(1.0);

    // These come from muscle.h
    //constructProperty_max_isometric_force(1000);
    //constructProperty_optimal_fiber_length(0.1);
    //constructProperty_tendon_slack_length(0.2);
    //constructProperty_max_contraction_velocity(10*optimal_fiber_length);
}

//_____________________________________________________________________________
// Set the data members of this muscle to their null values.
void VandenBogert2011Muscle::setNull()
{
    // no data members
    setAuthors("A. van den Bogert, B. Humphreys, C. Dembia");
    setReferences("van den Bogert A, Blana D, Heinrich D, Implicit methods for efficient musculoskeletal simulation and optimal control. Procedia IUTAM 2011; 2:297-316");

}


void VandenBogert2011Muscle::extendFinalizeFromProperties()
{
    Super::extendFinalizeFromProperties();

}

//--------------------------------------------------------------------------
// GET & SET Properties
//--------------------------------------------------------------------------
void VandenBogert2011Muscle::setFMaxTendonStrain(double fMaxTendonStrain) {
    set_tendon_strain_at_max_iso_force(fMaxTendonStrain); }
double VandenBogert2011Muscle::getFMaxTendonStrain() const {
     return get_tendon_strain_at_max_iso_force(); }

void VandenBogert2011Muscle::setFlWidth(double flWidth) {
     set_active_force_length_curve_width(flWidth); }
double VandenBogert2011Muscle::getFlWidth() const
    { return get_active_force_length_curve_width(); }

void VandenBogert2011Muscle::setFvAHill(double fvAHill) {
    set_force_velocity_hill_constant(fvAHill); }
double VandenBogert2011Muscle::getFvAHill() const {
     return get_force_velocity_hill_constant(); }

void VandenBogert2011Muscle::setFvmaxMultiplier(double fvMaxMultiplier) {
     set_force_velocity_max_lengthening_force_norm(fvMaxMultiplier); }
double VandenBogert2011Muscle::getFvmaxMultiplier() const {
     return get_force_velocity_max_lengthening_force_norm(); }

void VandenBogert2011Muscle::setDampingCoefficient(double dampingCoefficient) {
     set_fiber_damping(dampingCoefficient); }
double VandenBogert2011Muscle::getDampingCoefficient() const {
     return get_fiber_damping(); }

void VandenBogert2011Muscle::setNormFiberSlackLength(double
                                                     normFiberSlackLength) {
     set_fiber_slack_length_norm(normFiberSlackLength); }
double VandenBogert2011Muscle::getNormFiberSlackLength() const {
     return get_fiber_slack_length_norm(); }

void VandenBogert2011Muscle::setActivationTimeConstant(double activTimeConstant) {
     set_activation_time_constant(activTimeConstant); }
double VandenBogert2011Muscle::getActivationTimeConstant() const {
     return get_activation_time_constant(); }

void VandenBogert2011Muscle::setDeactivationTimeConstant(double
                                                    deactivationTimeConstant) {
     set_deactivation_time_constant(deactivationTimeConstant); }
double VandenBogert2011Muscle::getDeactivationTimeConstant() const {
     return get_deactivation_time_constant(); }

void VandenBogert2011Muscle::setPennAtOptFiberLength(double
                                                     pennAtOptFiberLength) {
    set_pennation_at_optimal_fiber_length(pennAtOptFiberLength); }
double VandenBogert2011Muscle::getPennAtOptFiberLength() const {
    return get_pennation_at_optimal_fiber_length(); }


void VandenBogert2011Muscle::setDefaultActivation(double
                                                     defaultActivation) {
    set_default_activation(defaultActivation); }
double VandenBogert2011Muscle::getDefaultActivation() const {
    return get_default_activation(); }

void VandenBogert2011Muscle::setDefaultFiberLength(double fiberLength)
{
    set_default_fiber_length(fiberLength);
}

double VandenBogert2011Muscle::getDefaultFiberLength() const
{   return get_default_fiber_length(); }


void VandenBogert2011Muscle::setProjFiberLengthNorm(SimTK::State& s, double projFibLenNorm)
            const
{
    //In other muscles this is setFiberLength
        setStateVariableValue(s, "projected_fiber_length_normalized",projFibLenNorm);
        markCacheVariableInvalid(s,"lengthInfo");
        markCacheVariableInvalid(s,"velInfo");
        markCacheVariableInvalid(s,"dynamicsInfo");
    }


// Set the activation
void VandenBogert2011Muscle::setActivation(SimTK::State& s,double activation)
const
//TODO: Add clamping?  Add ignore_activation_dynamics
{
    setStateVariableValue(s, "activation", activation);

    markCacheVariableInvalid(s,"velInfo");
    markCacheVariableInvalid(s,"dynamicsInfo");
}

// Get the Residual of muscle
// TODO:  This breaks naming convention in that it does not get an already
// evaluated constant, it actually performs the calculation.  Without using a
// state variable or adding to the cache variables, I think I have to do it this
// way?
SimTK::Vec2 VandenBogert2011Muscle::getResidual(const SimTK::State& s,
                                                double projFibVelNorm_guess,
                                                double activdot_guess,
                                                double excitation)
                                                const
{
    ImplicitResidual residualStruct = calcImplicitResidual(s,
                                                           projFibVelNorm_guess,
                                                           activdot_guess,
                                                           excitation, 0);

    SimTK::Vec2 residual = {residualStruct.forceResidual,
                            residualStruct.activResidual};
    return(residual);
}



//These are hacks because I am not using muscle.h mli cache variables, yet:
double  VandenBogert2011Muscle::getProjFiberLengthNorm(const SimTK::State& s) const
{
    return getStateVariableValue(s, "projected_fiber_length_normalized");
}
double  VandenBogert2011Muscle::getActivation(const SimTK::State& s) const
{
    return getStateVariableValue(s, "activation");
}



//==============================================================================
// Muscle.h Interface
//==============================================================================


double  VandenBogert2011Muscle::computeActuation(const SimTK::State& s) const
{

    const MuscleDynamicsInfo& mdi = getMuscleDynamicsInfo(s);

    setActuation(s, mdi.tendonForce);

    return mdi.tendonForce;
}



void VandenBogert2011Muscle::computeInitialFiberEquilibrium(SimTK::State& s)
const
{
    // TODO:  this is the same code as computeFiberEquilibriumAtZeroVelocity
    //Need to update for v~=0
    //Let's start with the default activation and the fiberVelocity=0

    computeFiberEquilibriumAtZeroVelocity(s);

}



void VandenBogert2011Muscle::
computeFiberEquilibriumAtZeroVelocity(SimTK::State& s) const
{

    double activation = getActivation(s);
    double projFibLenNorm = getProjFiberLengthNorm(s);
    SimTK::Vec2 lenghAndForce = calcFiberStaticEquilbirum(projFibLenNorm,
                                                        activation);

    setActuation(s, lenghAndForce[1]);
    setProjFiberLengthNorm(s,lenghAndForce[0]);

    //TODO:  Millard Muscle handles non-convergence here and zero muscle length.
    //     Need to consider adding.

}


void VandenBogert2011Muscle::
postScale(const SimTK::State& s, const ScaleSet& aScaleSet)
{

    GeometryPath& path = upd_GeometryPath();
    path.postScale(s, aScaleSet);

    if (path.getPreScaleLength(s) > 0.0) {
        double scaleFactor = getLength(s) / path.getPreScaleLength(s);
        upd_optimal_fiber_length() *= scaleFactor;
        upd_tendon_slack_length() *= scaleFactor;
        path.setPreScaleLength(s, 0.0);
    }
}

//==============================================================================
// MODELCOMPONENT INTERFACE REQUIREMENTS
//==============================================================================


void VandenBogert2011Muscle::extendConnectToModel(Model& model)
{
    Super::extendConnectToModel(model);
}


// Define new states and their derivatives in the underlying system
void VandenBogert2011Muscle::extendAddToSystem(SimTK::MultibodySystem& system)
const
{

    Super::extendAddToSystem(system);

    SimTK_ASSERT(isObjectUpToDateWithProperties(),
                 "VandenBogert2011Muscle: Muscle properties not up-to-date");

    addStateVariable("projected_fiber_length_normalized");
    addStateVariable("activation");

}

void VandenBogert2011Muscle::extendInitStateFromProperties(SimTK::State& s)
const
{
    Super::extendInitStateFromProperties(s);

    //Brad - Thelen does not do this here, but Millard does
    setActivation(s, getDefaultActivation());
    double projFiberLength = fiberLengthToProjectedLength(getDefaultFiberLength(), false);
    setProjFiberLengthNorm(s, projFiberLength/getOptimalFiberLength());
}

void VandenBogert2011Muscle::extendSetPropertiesFromState(const SimTK::State& s)
{
    Super::extendSetPropertiesFromState(s);

    //Brad - Thelen does not do this here, but Millard does
    setDefaultActivation(getStateVariableValue(s,"activation"));
    setDefaultFiberLength(getStateVariableValue(s,"projected_fiber_length_normalized"));

}


void VandenBogert2011Muscle::computeStateVariableDerivatives(const SimTK::State& s) const
{


    double excitation = getControl(s);

    // For now we will start with static guess
    // TODO: figure out how to use previous state as guess
    SimTK::Vector ydotInitialGuess(2);
    ydotInitialGuess[0] = 0.0;
    ydotInitialGuess[1] = 0.0;
    // Brad double projFibVelNorm = calcSolveMuscle(s, excitation, ydotInitialGuess);
    double projFibVelNorm = 0.05;


    setStateVariableDerivativeValue(s, "projected_fiber_length_normalized", projFibVelNorm );


    double adot =getActivationDerivative(s);

    setStateVariableDerivativeValue(s, "activation",  adot);

    cout << "excitation: " << excitation << endl;
    cout << "projFibVelNorm: " << projFibVelNorm << endl;
    cout << "adot: " << adot<< endl;

}


//=============================================================================
// COMPUTATION
//=============================================================================
class ImplicitSystemForwardEulerStep : public SimTK::OptimizerSystem {
public:

    //TODO:  Add back in option to use NR solver (and perform benchmark)
    /* Constructor class. Parameters passed are accessed in the objectiveFunc()
     * class. */
    ImplicitSystemForwardEulerStep(int numParameters, SimTK::State& s, double& u): SimTK::OptimizerSystem(numParameters), numControls(numParameters), si(s), excitation(u)
    {

        //VandenBogert2011Muscle::ImplicitResidual  Results1 = VandenBogert2011Muscle::calcImplicitResidual(s,0.0,0.0,0.0,0);

    }

    int objectiveFunc(const SimTK::Vector &new_yDotGuess, bool new_params,
                      SimTK::Real& f) const override {
        // No objective function
        f = 0;
        return 0;
    }

    int constraintFunc(const SimTK::Vector &new_yDotGuess, bool newParams,
                       SimTK::Vector& constraints)  const override {
        //The constraints function is the residuals

        double projFibVelNorm_guess = new_yDotGuess[0];
        double activdot = 0;  //BTH - big hack here.  Remove this

        VandenBogert2011Muscle::ImplicitResidual results = p_muscle->calcImplicitResidual(si,projFibVelNorm_guess,activdot,excitation,0);


        constraints[0] = results.forceResidual;
        //constraints[1] = results.activResidual;

        return 0;
    }

private:
    int numControls;
    SimTK::State& si;
    double excitation;
    SimTK::ReferencePtr<VandenBogert2011Muscle> p_muscle;
};



double VandenBogert2011Muscle::fiberLengthToProjectedLength
        (double fiberLength , bool argsAreNormalized) const {

    /* When argsAreNormalized= True the function expects normalized values and
     * returns normalized values.  The function though is internally structured
     * to utilize non-normalized values (therefore for it will convert
     * normalized arguments to non-normalized values, perform calculations,
     * and do then convert the return values to normalized
     * when argsAreNormalized= True.)*/

    double pennAtOptFiberLength = getPennAtOptFiberLength();
    double projFibLen = 0;


    // If there is pennation, compute fiberLength and muscleWidth from state s using the constant volume constraint: Lce*sin(p) = Lceopt*sin(popt)
    // Lce is dimensionless (normalized to Lceopt)
    // If pennation is zero, we can't do this because volume is zero, and fiberLength = projFibLen
    if (pennAtOptFiberLength < 0.01) return fiberLength;

    //else:
    double optimalFiberLength = getOptimalFiberLength();

    if (argsAreNormalized) fiberLength = fiberLength*optimalFiberLength;

    double muscleWidth = sin(pennAtOptFiberLength);
    //b is the fixed distance of the fiber perpindicular to the tendon (width)

    /*TODO: muscleWidth as cache variable
     * recalculated. Should be moved info muscleLengthInfo cache */

    // The fiber can not be shorter than b; if it is, the projected length
    // equation below will return a complex value.  It also physically can
    // not happen (the fiber being shorter than the fixed width)
    if (fiberLength >= muscleWidth) {
        projFibLen = sqrt(pow(fiberLength, 2) - pow(muscleWidth, 2));
    } else {
        projFibLen = SimTK::NaN;
    }  //TODO: Doing this for now, but need to
    // come up with a clamping scheme (Millard has one)

    if (argsAreNormalized) projFibLen = projFibLen/optimalFiberLength;


    return projFibLen;
};


//------------------------------------------------------------------------------
double VandenBogert2011Muscle::projFibLenToFiberLength (double projFibLen,
                                                        bool argsAreNormalized)
                                                        const {

    /* When argsAreNormalized= True the function expects normalized values and
     * returns normalized values.  The function though is internally structured
     * to utilize non-normalized values (therefore it will convert
     * normalized arguments to non-normalized values, perform calculations,
     * and do then convert the return values to normalized
     * when argsAreNormalized= True.)*/


    double pennAtOptFiberLength = getPennAtOptFiberLength();
    double fiberLength = 0;

    double optimalFiberLength = getOptimalFiberLength();

    if (argsAreNormalized) projFibLen = projFibLen *
                                        optimalFiberLength;  //Convert from Norm


    // If pennation is zero, we can't do this because volume is zero, and fiberLength = projFibLen
    if (pennAtOptFiberLength < 0.01) {
        fiberLength = projFibLen;
    } else {
        double static muscleWidth = sin(pennAtOptFiberLength);
        fiberLength = sqrt(pow(projFibLen, 2) + pow(muscleWidth, 2));
    }

    //Convert back to Norm if needed
    if (argsAreNormalized) fiberLength = fiberLength / optimalFiberLength;

    return fiberLength;
}

//------------------------------------------------------------------------------
double VandenBogert2011Muscle::fiberVelocityToProjFibVel
        (double fiberVelocity, double fiberLength, double projFiberLength,
         bool argsAreNormalized) const {

    /* When argsAreNormalized= True the function expects normalized values and
     * returns normalized values.  The function though is internally structured
     * to utilize non-normalized values (therefore for it will convert
     * normalized arguments to non-normalized values, perform calculations,
     * and do then convert the return values to normalized
     * when argsAreNormalized= True.)*/

    // Note that this normalized by fiber optimal length (not by max velocity)
    double optimalFiberLength = getOptimalFiberLength();

    if (argsAreNormalized) {
        fiberLength = fiberLength * getOptimalFiberLength();
        fiberVelocity = fiberVelocity * optimalFiberLength;
        projFiberLength = projFiberLength * optimalFiberLength;
    }

    double projFibVel = 0;

    if (fiberVelocity != 0) {
        projFibVel = fiberVelocity / cos(projFiberLength / fiberLength);

    } else {
        projFibVel = 0; }

    if (argsAreNormalized) projFibVel = projFibVel / optimalFiberLength;

    return projFibVel;
};

//------------------------------------------------------------------------------
double VandenBogert2011Muscle::fiberVelocityToProjFibVel
        (double fiberVelocity, double fiberLength, bool argsAreNormalized)
                const {

    //overload method when projFiberLength is not known/provided

    double projFiberLength = fiberLengthToProjectedLength(fiberLength,
                                                          argsAreNormalized);
    double projFiberVelocity = fiberVelocityToProjFibVel (fiberVelocity,
                                                          fiberLength,
                                                          projFiberLength,
                                                          argsAreNormalized);
    return projFiberVelocity;
};


//------------------------------------------------------------------------------
double VandenBogert2011Muscle::projFibVelToFiberVelocity(double projFibVel,
                                                         double fiberLength,
                                                         double projFiberLength,
                                                         bool argsAreNormalized)
                                                         const {

    /* When argsAreNormalized= True the function expects normalized values and
 * returns normalized values.  The function though is internally structured
 * to utilize non-normalized values (therefore for it will convert
 * normalized arguments to non-normalized values, perform calculations,
 * and do then convert the return values to normalized
 * when argsAreNormalized= True.)*/

// Note that this normalized by fiber optimal length (not by max velocity)
    double optimalFiberLength = getOptimalFiberLength();

    if (argsAreNormalized) {
        fiberLength = fiberLength * getOptimalFiberLength();
        projFibVel = projFibVel * optimalFiberLength;
        projFiberLength = projFiberLength * optimalFiberLength;
    }

    double fiberVelocity = 0;

    if (projFibVel != 0) {
        fiberVelocity = projFibVel * cos(projFiberLength / fiberLength);
    } else {
        fiberVelocity = 0; }

    if (argsAreNormalized) fiberVelocity = fiberVelocity / optimalFiberLength;

    return fiberVelocity;}

//------------------------------------------------------------------------------
double VandenBogert2011Muscle::projFibVelToFiberVelocity
        (double projFibVel, double projFibLength, bool argsAreNormalized)
                        const {

    //overload method when fiberLength is not known/provided

    double fiberLength = projFibLenToFiberLength(fiberLength, argsAreNormalized);
    double fiberVelocity = projFibVelToFiberVelocity(projFibVel, fiberLength,
                                                     projFibLength,
                                                     argsAreNormalized);

    return fiberVelocity;
};



//------------------------------------------------------------------------------
VandenBogert2011Muscle::ImplicitResidual
VandenBogert2011Muscle::calcImplicitResidual(const SimTK::State& s,
                                             double projFibVelNorm_guess,
                                             double activdot_guess,
                                             double excitation,
                                             bool returnJacobians)
const
{
    // Overload method for state as parameters

    VandenBogert2011Muscle::ImplicitResidual results = calcImplicitResidual(
            getLength(s), getProjFiberLengthNorm(s), getActivation(s), projFibVelNorm_guess,
            activdot_guess, excitation, returnJacobians);

    return results;
}


//------------------------------------------------------------------------------
VandenBogert2011Muscle::ImplicitResidual
VandenBogert2011Muscle::calcImplicitResidual(SimTK::Vec2 y,
                                             SimTK::Vec2 ydot_guess,
                                             double muscleLength,
                                             double excitation, bool returnJacobians)
const
{
    // Overload method for state vectors as parameters
    VandenBogert2011Muscle::ImplicitResidual results = calcImplicitResidual(
            muscleLength, y[0], y[1], ydot_guess[0], ydot_guess[1], excitation,
            returnJacobians);
    return results;

}

//------------------------------------------------------------------------------
VandenBogert2011Muscle::ImplicitResidual VandenBogert2011Muscle::
calcImplicitResidual(double muscleLength, double projFibLenNorm, double activ,
                     double projFibVelNorm, double activdot, double excitation,
                     bool returnJacobians) const {


    //TODO: Match symbols to doxygen diagram and Add symbolic equations comments
    //TODO: May want to make this into a separate function

    // -------------------------Parameters----------------------------//
    //F_{iso} Maximum isometric force that the fibers can generate
    double maxIsoForce = getMaxIsometricForce();

    //u_{max} (dimensionless) strain in the series elastic element at load of
    //      maxIsometricForce
    double fMaxTendonStrain = getFMaxTendonStrain(); //Strain in the series
    //      elastic element at load of Fmax

    //W (dimensionless) width parameter of the force-length relationship of
    //      the muscle fiber
    double fl_width = getFlWidth();

    //AHill (dimensionless) Hill parameter of the force-velocity relationship
    double fv_AHill = getFvAHill();   //Hill parameter of the f-v relationship

    //FV_{max} (dimensionless) maximal eccentric force
    double fv_maxMultiplier = getFvmaxMultiplier();
    //Maximal eccentric force multiplier

    // L_{opt}  (m) Optimal Length of Contractile Element;
    double optFiberLength = getOptimalFiberLength();
    //Optimal Length of Contractile Element

    //b (s/m) damping coefficient of damper parallel to the fiber
    // (normalized to maxIsometricForce)
    double dampingCoeff = getDampingCoefficient();

    //L_{slack,fiber} Slack length of the parallel elastic element, divided by
    //          Lceopt
    double fiberSlackLengthNorm = getNormFiberSlackLength();

    //L_{slack,tendon} (m) slack length of the tendon
    double tendonSlackLength = getTendonSlackLength();

    //T_{act} (s) Activation time
    double activTimeConstant = getActivationTimeConstant();

    //T_{deact} (s) Deactivation time
    double deactivationTimeConstant = getDeactivationTimeConstant();

    //phi_{opt} pennation at Optimal Fiber Length
    double pennAtOptFiberLength = getPennAtOptFiberLength();


    // constants derived from the muscle parameters


    // Jacobian Matrices
    SimTK::Mat22 df_dy;
    SimTK::Mat22 df_dydot;


    //-------Convert projFiberLength & Velocity to fiberLength & Velocity------/
    double fiberLengthNorm;
    double dfiberLength_dprojFibLen = SimTK::NaN;
    double dcosPenn_dprojFibLen = SimTK::NaN;

    fiberLengthNorm=projFibLenToFiberLength(projFibLenNorm,true);

    SimTK::Vec3 cosPennAndDeriv=penMdl_calcCosPennationAngle(projFibLenNorm,fiberLengthNorm,returnJacobians);
    double cosPenn = cosPennAndDeriv[0];

    if (returnJacobians) {
        dfiberLength_dprojFibLen = cosPennAndDeriv[1];
        dcosPenn_dprojFibLen = cosPennAndDeriv[2];
    }


    // Compute fiberVelocity and its derivatives wrt projFibLen and projFibVel
    double fiberLengtheningVelocityNorm = projFibVelNorm*cosPenn;
    // This looks incorrect on first inspection (that it should be /cosPenn)
    // If you take the time derivative using the chain rule of:
    //      l_m=(l_p^2+h^2)^0.5
    // you get l_mdot=l_pdot * l_p/l_m = l_pdot * cos(penn)


    //--- F1 is the normalized isometric force-length relationship at maximum
    //                                               activation--------------//
    SimTK::Vec2 fiberForceLengthMultiplierAndDeriv = calcFiberActiveForceLengthMultiplier(fiberLengthNorm, returnJacobians);
    double F1 = fiberForceLengthMultiplierAndDeriv[0];
    double dF1_dfiberLengthNorm = fiberForceLengthMultiplierAndDeriv[0];


    //-------- F2 is the dimensionless force-velocity relationship -------//
    SimTK::Vec8 fv = calcFiberActiveForceVelocityMultiplier(fiberLengtheningVelocityNorm, cosPenn, fiberLengthNorm, returnJacobians);

    double F2 = fv[0];
    double dF2_dactiv =SimTK::NaN;
    double dF2_dprojFibLenNorm = SimTK::NaN;
    double dF2_dprojFibVelNorm = SimTK::NaN;
    if (returnJacobians) {
        double dF2_dfiberLengtheningVelocityNorm = fv[1];
        dF2_dactiv = fv[2];
        dF2_dprojFibVelNorm = fv[3];
        dF2_dprojFibLenNorm = fv[4];
        double dfiberLengtheningVelocityNorm_dprojFibVelNorm = fv[6];
        double dfiberLengtheningVelocityNorm_dprojFibLenNorm = fv[7];
    }

    //------F3 is the dimensionless fiber (PEE) force (in units of Fmax)------//
    SimTK::Vec2 fiberPassiveForceLengthMultiplierAndDeriv = calcFiberPassiveForceLengthMultiplier(fiberLengthNorm, returnJacobians);
    double F3 = fiberPassiveForceLengthMultiplierAndDeriv[0];
    double dF3_dfiberLengthNorm = SimTK::NaN;
    if (returnJacobians) {
        dF3_dfiberLengthNorm = fiberPassiveForceLengthMultiplierAndDeriv[1];
    }

    //--------F4 is the dimensionless SEE force (in units of Fmax)----------//
    SimTK::Vec4 tendonForceMultiplierAndDeriv = calcTendonForceMultiplier(muscleLength, projFibLenNorm, returnJacobians);
    double F4 = tendonForceMultiplierAndDeriv[0];
    double dF4_dprojFibLenNorm = SimTK::NaN;
    double dF4_dmuscleLength = SimTK::NaN;
    if (returnJacobians) {
        dF4_dprojFibLenNorm = tendonForceMultiplierAndDeriv[1];
        dF4_dmuscleLength = tendonForceMultiplierAndDeriv[2];
    }

    //-- F5 is viscous damping parallel to the CE (0.001 of Fmax at 1 Lceopt/s)
    // to  ensure that df/dLcedot is never zero-----------//
    double F5 = dampingCoeff * projFibVelNorm ;
    double dF5_dprojFibVelNorm  = dampingCoeff;

    // ---------Calculate the Muscle Force Residual ---------------------- //
    //The muscle dynamics equation: f = Fsee - (a*Fce - Fpee)*cos(Penn) -
    //                                                          Fdamping = 0
    //   (Reminder the units of fRes are (N/N), needs to be *Fmax to get to N)
    double fRes = F4 - (activ * F1 * F2 + F3)*cosPenn - F5;

if (getDebugLevel()>0){
    cout << "-------------------------" << endl;
    cout << "activ: " << activ << endl;
    cout << "F1 (FL): " << F1 << endl;
    cout << "F2 (FV): " << F2 << endl;
    cout << "F3 (PEE) : " << F3 << endl;
    cout << "F4 (SEE): " << F4 << endl;
    cout << "F5 (Damping): " << F5 << endl;
    cout << "fRes: " << fRes << endl;
    cout << "------------------" << endl;
    cout << "cosPenn: " << cosPenn << endl;
    cout << "muscleLength: " << muscleLength << endl;
    cout << " SEE" << endl;
    //cout << "   kSEE: " << kSEE << endl;
    //cout << "   kSEE2: " << kSEE2 << endl;
    //cout << "   elongationTendon: " << elongationTendon << endl;
    cout << "   tendonSlackLength: " << tendonSlackLength << endl;
    cout << " PEE" << endl;
    cout << "   fiberLengthNorm: " << fiberLengthNorm << endl;
    cout << "   optFiberLength: " << optFiberLength << endl;
    cout << "   projFibLenNorm: " << projFibLenNorm << endl;
    //cout << "   kPEENorm: " << kPEENorm << endl;
    //cout << "   kPEE2Norm: " << kPEE2Norm << endl;
    //cout << "   elongationFiberNorm: " <<  elongationFiberNorm << endl;
    cout << "" << endl;
}

    // --------------- Force in tendon (SEE) is maxIsoForce*F4 -------------- //
    double Fsee = maxIsoForce * F4;

    // TODO: I don't think we need to implement dFsee_.....
    /*if (nargout > 1)
        dFsee_dLm  = Fmax * dF4_dLm;
    dFsee_dLce = Fmax * dF4_dLce;
    dFsee_dy   = [0;0;(-d/L)*dFsee_dLm; 0;dFsee_dLce;0;0];
     % Fsee is a function of y(3) & y(5)
    end*/

    //----------------------Activation dynamics equation-------------------//

    double activationResidual = activdot - calcActivationDerivative(activ,excitation);
    SimTK::Vec2 df_du;
    double dActRes_dactiv = 0;
    double dActRes_dactivdot = 0;

    if (returnJacobians) {
        dActRes_dactiv = (excitation / activTimeConstant + (1 - excitation) / deactivationTimeConstant);
        dActRes_dactivdot = 1;

        df_du[0] = 0;
        df_du[1] = -(excitation / activTimeConstant + (1 - excitation) / deactivationTimeConstant )
                       - (excitation - activ) * (1 / activTimeConstant - 1 / deactivationTimeConstant );
    }

    //---------------------Assemble Jacobians---------------------------//
    if (returnJacobians) {
        double dfRes_dactiv = -(F1*F2 + activ*F1*dF2_dactiv )*cosPenn;


        double dF3_dprojFibLenNorm = dF3_dfiberLengthNorm  * dfiberLength_dprojFibLen;

        double dfRes_dprojFibLengthNorm = dF4_dprojFibLenNorm -
              (activ*(dF1_dfiberLengthNorm*F2 + F1*dF2_dprojFibLenNorm) +
               dF3_dprojFibLenNorm) * cosPenn - (activ*F1*F2 + F3) *
                dcosPenn_dprojFibLen;

        double dfRes_dprojFibVelNorm =
                - activ*F1*dF2_dprojFibVelNorm - dF5_dprojFibVelNorm;

        //y=(projFiberLength,activation)  <-States

        //df_fy is 2x2:  Where columns are states:
        //          [projFiberLengthNorm, activation]
        //              Rows are:
        //          [force residual; Activation Residual]

        // Row 1 - df_dy  (force)
        df_dy[0][0] = dfRes_dprojFibLengthNorm;
        df_dy[0][1] = dfRes_dactiv;
        // Row 2 - df_dy (activation)
        df_dy[1][0] = 0;
        df_dy[1][1] = dActRes_dactiv ;


        // Row 1 - df_dydot  (force)
        df_dydot[0][0] = dfRes_dprojFibVelNorm;
        df_dydot[0][1] = 0;
        // Row 2 - df_dydot (activation)
        df_dydot[1][0] = 0;
        df_dydot[1][1] = dActRes_dactivdot;

    }



    VandenBogert2011Muscle::ImplicitResidual results;

    results.forceResidual = fRes;
    results.activResidual = activationResidual;
    results.forceTendon = Fsee;
    results.df_dy = df_dy;
    results.df_dydot = df_dydot;
    results.df_du = df_du;
    results.df_dmuscleLength = dF4_dmuscleLength;
    results.F1 = F1;    //Output force components for troubleshooting
    results.F2 = F2;
    results.F3 = F3;
    results.F4 = F4;
    results.F5 = F5;

return results; }



double VandenBogert2011Muscle::
calcActivationDerivative(double activation, double excitation) const
{

    double activationDerivative=(excitation - activation) *
    (excitation / getActivationTimeConstant() + (1 - excitation) / getDeactivationTimeConstant());

    return activationDerivative;
}


double VandenBogert2011Muscle::
getActivationDerivative(const SimTK::State& s) const
{
    double activationDerivative =
            calcActivationDerivative(getActivation(s),getExcitation(s));

    return activationDerivative;
}





//==============================================================================
// MUSCLE INTERFACE REQUIREMENTS -- MUSCLE LENGTH INFO
//==============================================================================
void VandenBogert2011Muscle::calcMuscleLengthInfo(const SimTK::State& s,
                                                        MuscleLengthInfo& mli) const
{
    // Get musculotendon actuator properties.
    //double maxIsoForce    = getMaxIsometricForce();
    double optFiberLength = getOptimalFiberLength();
    double tendonSlackLen = getTendonSlackLength();

    try {

        double projFibLengthNorm = getStateVariableValue(s, "projected_fiber_length_normalized");


        mli.normFiberLength = projFibLenToFiberLength(projFibLengthNorm,true);
        mli.fiberLength = mli.normFiberLength * optFiberLength ;
        cout << "projFibLengthNorm : " << projFibLengthNorm  << endl;
        cout << "mli.fiberLength : " << mli.fiberLength  << endl;
        cout << "mli.normFiberLength: " << mli.normFiberLength << endl;
        cout << "optFiberLength: " << optFiberLength << endl;


        SimTK::Vec3 cosPennAndDeriv = penMdl_calcCosPennationAngle(projFibLengthNorm, mli.normFiberLength, true);

        mli.cosPennationAngle = cosPennAndDeriv[0];
        mli.pennationAngle = acos(mli.cosPennationAngle);

        //TODO: Do we want to report the derivatives?
        mli.userDefinedLengthExtras.resize(2);
        mli.userDefinedLengthExtras[0] = cosPennAndDeriv[1]; //dfiberLength_dprojFibLen
        mli.userDefinedLengthExtras[1] = cosPennAndDeriv[2]; //dcosPenn_dprojFibLen

        mli.sinPennationAngle = sin(mli.pennationAngle);
        mli.fiberLengthAlongTendon = mli.fiberLength * mli.cosPennationAngle;

        mli.tendonLength      = penMdl_calcTendonLength(projFibLengthNorm, getLength(s));

        mli.normTendonLength  = mli.tendonLength / tendonSlackLen;
        mli.tendonStrain      = mli.normTendonLength - 1.0;

        //TODO: Do we want to report the derivatives?
        SimTK::Vec2 fiberPassiveForceLengthMultiplierAndDeriv = calcFiberPassiveForceLengthMultiplier(mli.normFiberLength, false);
        mli.fiberPassiveForceLengthMultiplier = fiberPassiveForceLengthMultiplierAndDeriv[0];

        SimTK::Vec2 fiberActiveForceLengthMultiplierAndDeriv = calcFiberActiveForceLengthMultiplier(mli.normFiberLength, false);
        mli.fiberActiveForceLengthMultiplier = fiberActiveForceLengthMultiplierAndDeriv[0];


    } catch(const std::exception &x) {
        std::string msg = "Exception caught in VandenBogert2011Muscle::"
                                  "calcMuscleLengthInfo from " + getName() + "\n"
                          + x.what();
        throw OpenSim::Exception(msg);
    }
}




//==============================================================================
// MUSCLE INTERFACE REQUIREMENTS -- MUSCLE POTENTIAL ENERGY INFO
//==============================================================================
void VandenBogert2011Muscle::
        calcMusclePotentialEnergyInfo(const SimTK::State& s,
                                      MusclePotentialEnergyInfo& mpei) const
{
try {
    // Get the quantities that we've already computed.
    const MuscleLengthInfo &mli = getMuscleLengthInfo(s);

    double projFibLengthNorm = getProjFiberLengthNorm(s);
    mpei.fiberPotentialEnergy = calcFiberPotentialEnergy(projFibLengthNorm);
    mpei.tendonPotentialEnergy = calcTendonPotentialEnergy(projFibLengthNorm,
                                                           mli.fiberLength);


    mpei.musclePotentialEnergy = mpei.fiberPotentialEnergy +
                                 mpei.tendonPotentialEnergy;

    } catch(const std::exception &x) {
        std::string msg = "Exception caught in VandenBogert2011Muscle::"
                                  "calcMusclePotentialEnergyInfo from " + getName() + "\n"
                          + x.what();
        throw OpenSim::Exception(msg);
    }
}





//==============================================================================
// MUSCLE INTERFACE REQUIREMENTS -- FIBER VELOCITY INFO
//==============================================================================
void VandenBogert2011Muscle::calcFiberVelocityInfo(const SimTK::State& s, FiberVelocityInfo& fvi) const

{
    try {


        // Get the quantities that we've already computed.
        const MuscleLengthInfo &mli = getMuscleLengthInfo(s);

        // Get the static properties of this muscle.
        //double dlenMcl   = getLengtheningSpeed(s);
        //double optFibLen = getOptimalFiberLength();


        double projFiberVelocityNorm = getProjFiberVelNorm(s);

        double optimalFiberLength = getOptimalFiberLength();

        double cosPenn = mli.cosPennationAngle;


        fvi.fiberVelocityAlongTendon     = projFiberVelocityNorm * optimalFiberLength;
        fvi.fiberVelocity = fvi.fiberVelocityAlongTendon*cosPenn;
        // This should be *cosPenn; see notes in calcImplicitResidual()
        fvi.normFiberVelocity = fvi.fiberVelocity / optimalFiberLength;


        fvi.pennationAngularVelocity     = penMdl_calcPennationAngularVelocity(
                mli.cosPennationAngle, mli.sinPennationAngle,
                mli.fiberLength, fvi.fiberVelocityAlongTendon);

        //TODO:  Verify this as Millard/MuscleFixedWidthPennationModel seems to
        // be much more complicated than this.
        fvi.tendonVelocity =  getLengtheningSpeed(s) - fvi.fiberVelocity;
        fvi.normTendonVelocity = fvi.tendonVelocity / getTendonSlackLength();


        SimTK::Vec8 fv = calcFiberActiveForceVelocityMultiplier(fvi.normFiberVelocity,
                                                 mli.cosPennationAngle,
                                                 mli.normFiberLength, false);

        fvi.fiberForceVelocityMultiplier = fv[0];

    } catch(const std::exception &x) {
        std::string msg = "Exception caught in VandenBogert2011Muscle::"
                                  "calcFiberVelocityInfo from " + getName() + "\n"
                          + x.what();
        throw OpenSim::Exception(msg);
    }
}

//==============================================================================
// MUSCLE INTERFACE REQUIREMENTS -- MUSCLE DYNAMICS INFO
//==============================================================================
void VandenBogert2011Muscle::
calcMuscleDynamicsInfo(const SimTK::State& s, MuscleDynamicsInfo& mdi) const {

    try {
        // Get the quantities that we've already computed.
        const MuscleLengthInfo &mli = getMuscleLengthInfo(s);
        const FiberVelocityInfo &mvi = getFiberVelocityInfo(s);

        // Get the properties of this muscle.
        double tendonSlackLen = getTendonSlackLength();
        double optFiberLen = getOptimalFiberLength();
        double maxIsometricForce = getMaxIsometricForce();

        // Compute the stiffness of the muscle fiber.
        SimTK_ERRCHK_ALWAYS(mli.fiberLength > SimTK::SignificantReal,
                            "calcMuscleDynamicsInfo",
                            "The muscle fiber has a length of 0, causing a singularity");


        double activation = getControl(s);


        SimTK::Vec8 fiberActiveForceVelocityMultiplier =
                calcFiberActiveForceVelocityMultiplier(mvi.normFiberVelocity,
                                                       mli.cosPennationAngle,
                                                       mli.normFiberLength,
                                                       false);

        SimTK::Vec2 fiberForceLengthMultiplier =
                calcFiberActiveForceLengthMultiplier(mli.normFiberLength,
                                                     true);

        SimTK::Vec2 fiberPassiveForceLengthMultiplier =
                calcFiberPassiveForceLengthMultiplier(mli.normFiberLength,
                                                      false);

        SimTK::Vec2 fiberDampingMultiplier =
                calcFiberDampingMultiplier(
                        mvi.fiberVelocityAlongTendon/optFiberLen, false);

        SimTK::Vec4 tendonForceMultiplier = calcTendonForceMultiplier(
                getLength(s), mli.fiberLengthAlongTendon/optFiberLen,
                true);


        double fibActFvMult = fiberActiveForceVelocityMultiplier[0];
        double fibActFlMult = fiberForceLengthMultiplier[0];
        double fibPassFLMult = fiberPassiveForceLengthMultiplier[0];
        double fibDampMult = fiberDampingMultiplier[0];
        double tendFMult = tendonForceMultiplier[0];


        double activeFiberForceNorm =
                activation * fibActFvMult * fibActFlMult;
        double fiberForceAlongTendonNorm =
                (activeFiberForceNorm + fibPassFLMult) *
                mli.cosPennationAngle - fibDampMult;

        mdi.activation = activation;
        mdi.fiberForceAlongTendon =
                fiberForceAlongTendonNorm * maxIsometricForce;

        //TODO: This will cause a singularity at cosPenn=90, but this is only in the mdi reporting
        //Need to think of how to handle fiberForce
        mdi.normFiberForce =
                fiberForceAlongTendonNorm / mli.cosPennationAngle;
        mdi.fiberForce = mdi.normFiberForce * maxIsometricForce;

        mdi.activeFiberForce = activeFiberForceNorm * maxIsometricForce;
        mdi.passiveFiberForce = fibPassFLMult * maxIsometricForce;
        mdi.tendonForce = tendFMult * maxIsometricForce;
        mdi.normTendonForce = tendFMult;


        double dfibActFlMult_dfiberLengthNorm = fiberForceLengthMultiplier[1];
        //double dfibPassFLMult_dfiberLengthNorm = fiberPassiveForceLengthMultiplier[1];

        mdi.fiberStiffness =
                (activation * dfibActFlMult_dfiberLengthNorm *
                 fibActFvMult + dfibActFlMult_dfiberLengthNorm) *
                (maxIsometricForce / optFiberLen);

        // $\frac{dF_{mAT}}{dl_{ceAT}} =F_{m} \frac{sin^2(penn)}{l_{ce}} +  \frac{dF_{m}}{dl_{ce}}  {cos^2(penn)}$
        mdi.fiberStiffnessAlongTendon =
                mdi.activeFiberForce * pow(mli.sinPennationAngle, 2) /
                mli.fiberLength +
                mdi.fiberStiffnessAlongTendon *
                pow(mli.cosPennationAngle, 2);

        mdi.tendonStiffness = tendonForceMultiplier[4] * getMaxIsometricForce();

        //TODO:  this seems odd.  From Millard2012EquilbriumMuscle:
        // Ke, muscleStffness = (dFmAT_dlceAT*dFt_dtl)/(dFmAT_dlceAT+dFt_dtl);
        // Why is the AT used for the fiber.  Since muscle length is along the
        // fiber, should it be along the fiber?  Copying what Millard:
        // mdi.muscleStiffness = (mdi.fiberStiffnessAlongTendon*mdi.tendonStiffness) / (mdi.fiberStiffnessAlongTendon + mdi.tendonStiffness);
        // but I think this should be:
        mdi.muscleStiffness = (mdi.fiberStiffness*mdi.tendonStiffness) / (mdi.fiberStiffness + mdi.tendonStiffness);

        // TODO:  Verify power entires and clean up
        //Adding these so I can copy in place the code from Millard Power Entries.  I have not verified the Millard power equations.
        double fiso = getMaxIsometricForce();
        double p1Fm = fiso * fibActFlMult;
        double p2Fm = fiso * fibDampMult;
        double fse = mdi.normTendonForce;

        double dFibPEdt     = p1Fm*mvi.fiberVelocity; //only conservative part
        //of passive fiber force
        double dTdnPEdt     = fse*fiso*mvi.tendonVelocity;
        double dFibWdt      = -(mdi.activeFiberForce+p2Fm)*mvi.fiberVelocity;
        double dmcldt       = getLengtheningSpeed(s);
        double dBoundaryWdt = mdi.tendonForce*dmcldt;

        // Populate the power entries.
        mdi.fiberActivePower  = dFibWdt;
        mdi.fiberPassivePower = -(dFibPEdt);
        mdi.tendonPower       = -dTdnPEdt;
        mdi.musclePower       = -dBoundaryWdt;


    } catch(const std::exception &x) {
    std::string msg = "Exception caught in VandenBogert2011Muscle::"
                              "calcMuscleDynamicsInfo from " + getName() + "\n"
                      + x.what();
    cerr << msg << endl;
    throw OpenSim::Exception(msg);
    }

}



double VandenBogert2011Muscle::getProjFiberVelNorm(const SimTK::State& s) const
{

    //TODO:  - Should make this use the current state derivative value if
    // possible.  Does it need to recalcuate?  It is called by the velocity
    // cache update. But I guess it would solve quickly if the muscle is
    // balanced.

    double excitation = getExcitation(s);
    SimTK::Vector projFibVelNormGuess(1);

    projFibVelNormGuess[0]=0;
    double projFiberVelocityNorm = calcSolveMuscle(s, excitation, projFibVelNormGuess);

    return projFiberVelocityNorm;
}

//------------------------------------------------------------------------------
VandenBogert2011Muscle::ImplicitResidual
VandenBogert2011Muscle::calcJacobianByFiniteDiff(double muscleLength,
                                                 double projFibLenNorm,
                                                 double activ,
                                                 double projFibVelNorm,
                                                 double activdot, double excitation,
                                                 double stepSize) const
{
    SimTK::Vec2 y;
    SimTK::Vec2 ydot;
    y[0]=projFibLenNorm;
    y[1]=activ;
    ydot[0]=projFibVelNorm;
    ydot[1]=activdot;

    // Overload method for state vectors as parameters
    VandenBogert2011Muscle::ImplicitResidual results =
            calcJacobianByFiniteDiff(y, ydot, muscleLength, excitation, stepSize );

    return results;
}


//------------------------------------------------------------------------------
VandenBogert2011Muscle::ImplicitResidual    VandenBogert2011Muscle::
calcJacobianByFiniteDiff(SimTK::Vec2 y, SimTK::Vec2 ydot, double muscleLength,
                         double excitation, double stepSize ) const {



    //TODO: In the optimizer class, you do not have to supply the Jacobian
    // function.  The must mean there is code in Simbody that does finite diff.
    // Should probably look at calling that in place of this code.

    // Jacobian Matrices
    SimTK::Mat22 df_dy;
    SimTK::Mat22 df_dydot;
    SimTK::Vec2 df_du;

    VandenBogert2011Muscle::ImplicitResidual opPoint;
    opPoint = calcImplicitResidual (y,ydot,muscleLength,excitation,0);
    double opForceResidual = opPoint.forceResidual;
    double opActivResidual = opPoint.activResidual;

    VandenBogert2011Muscle::ImplicitResidual del;

    //----------df_dy------------//
    SimTK::Vec2 dh = {stepSize,0};
    del = calcImplicitResidual(y+dh,ydot,muscleLength,excitation,0);
    df_dy[0][0] = (del.forceResidual-opForceResidual)/stepSize;
    df_dy[1][0] = (del.activResidual-opActivResidual)/stepSize;

    dh = {0,stepSize};
    del = calcImplicitResidual(y+dh,ydot,muscleLength,excitation,0);
    df_dy[1][0] = (del.forceResidual-opForceResidual)/stepSize;
    df_dy[1][1] = (del.activResidual-opActivResidual)/stepSize;

    //----------df_dydot------------//
    dh = {stepSize,0};
    del = calcImplicitResidual(y,ydot+dh,muscleLength,excitation,0);
    df_dydot[0][0] = (del.forceResidual-opForceResidual)/stepSize;
    df_dydot[1][0] = (del.activResidual-opActivResidual)/stepSize;

    dh = {0,stepSize};
    del = calcImplicitResidual(y,ydot+dh,muscleLength,excitation,0);
    df_dydot[1][0] = (del.forceResidual-opForceResidual)/stepSize;
    df_dydot[1][1] = (del.activResidual-opActivResidual)/stepSize;

    //----------df_du----------------//
    del = calcImplicitResidual(y,ydot,muscleLength,excitation+stepSize,0);
    df_du[0] = (del.forceResidual-opForceResidual)/stepSize;
    df_du[1] = (del.activResidual-opActivResidual)/stepSize;

/*
    for (int i =0; i<=1; i=i+1) {
        yTemp[i]=y[i]+h;
        del = calcImplicitResidual(yTemp,ydot,u,0);
        df_dy[0][i]=(del.forceResidual-opForceResidual)/h; //force residual
        df_dy[1][i]=(del.activResidual-opActivResidual)/h; //activ residual


        ydotTemp[i]=y[i]+h;
        del = calcImplicitResidual(y,ydotTemp,u,0);

        df_dydot[0][i]=(del.forceResidual-opForceResidual)/h; //force residual
        df_dydot[1][i]=(del.activResidual-opActivResidual)/h; //activ residual

        yTemp=y;
        ydotTemp=ydot;

    }*/


    //del = calcImplicitResidual(y,ydot,u+h,0);
    //double df_du = (del.forceResidual-opForceResidual)/h;

    opPoint.df_dy = df_dy;
    opPoint.df_dydot = df_dydot;
    opPoint.df_du = df_du;

    return opPoint;
}


//------------------------------------------------------------------------------
// Just a quick convenience method to get the force residual
// under static conditions
SimTK::Vec3 VandenBogert2011Muscle::calcFiberStaticEquilibResidual(
        double projFibLenNorm, double muscleLength, double activ) const
{
    VandenBogert2011Muscle::ImplicitResidual results = calcImplicitResidual(
            muscleLength, projFibLenNorm, activ, 0.0, 0.0, activ, 1.0);

    SimTK::Vec3 resAndDerivative;

    resAndDerivative[0] = results.forceResidual;
    resAndDerivative[1] = results.df_dy[0][0];
    resAndDerivative[2] = results.forceTendon;

    return resAndDerivative;
}

//------------------------------------------------------------------------------
SimTK::Vec2 VandenBogert2011Muscle::calcFiberStaticEquilbirum(
        double muscleLength, double activ) const {

    // This code assumes that the muscle lengthing speed is 0.
    // It utilizes a Newton-Raphson Method to root solve to obtain the
    //fiber length.


    //TODO: Code is not optimized.  Specifically the number of calls to
    //calcImplicitResidual can be reduced

    //TODO: Generalize with a Lambda function (will need help with that).
    //TODO: calcImplicitResidual really only needs to calculate df_ds
    //       (single element) for this function


    double tol = 1e-8; //TODO : Look at changing to proportional of eps
    double a = 0;
    double b = 10;  //10 muscle lengths (more will slow convergence by adding steps)

    double x = (a + b) / 2.0;
    double dx = 2 * tol;

    int neval = 0;

    SimTK::Vec3 forceResAndDerivative;

    //Perform a Newton Line Search
    while ((abs(dx) >= tol) && (neval < 100)) {

        neval++;

        // Set a to be lower value and b to be upper value
        a = min(a, b);
        b = max(a, b);


        forceResAndDerivative =
                calcFiberStaticEquilibResidual(x, muscleLength, activ);
        double fx = forceResAndDerivative[0];

        //After the 1st iteration, use the new guess as a new upper or lower bound
        if (neval > 1) {
            forceResAndDerivative =
                    calcFiberStaticEquilibResidual(a, muscleLength, activ);
            double funcA = forceResAndDerivative[0];

            if ((funcA * fx) > 0) {
                a = x;
            } else {
                b = x;
            }

            forceResAndDerivative =
                    calcFiberStaticEquilibResidual(x, muscleLength, activ);
            double dfx = forceResAndDerivative[1];
            //double forceRes = forceResAndDerivative[0];

            dx = -fx / dfx;
            double xdx = x - dx;

            bool inInterval = ((xdx >= a) && (xdx <= b));

            forceResAndDerivative =
                    calcFiberStaticEquilibResidual(xdx, muscleLength, activ);
            bool largeDeriv = abs(forceResAndDerivative[0]) > (0.5 * abs(fx));

            if (~inInterval || largeDeriv) {
                x = (a + b) / 2;
                dx = (a - b) / 2;
            } else {
                x = xdx;
            }

            // TODO:  Need to handle condition when number of loop iterations reaches
            // neval limit. See Millard2012 Muscle Error Handling
        }
    }
    SimTK::Vec2 vout;
    vout[0] = x;   //projFiberLengthNorm
    vout[1] = forceResAndDerivative[2];  // muscleForce
    return vout;
};


//------------------------------------------------------------------------------
//Calculate the ydot values to drive the residuals to 0 and "balance" the muscle

double VandenBogert2011Muscle::calcSolveMuscle(const SimTK::State& s,
          double excitation, SimTK::Vector projFibVelNormGuess) const {

//SimTK::Vector VandenBogert2011Muscle::calcSolveMuscle(const SimTK::State& s,
//         double activ, SimTK::Vector yDotInitialGuess)  {

    SimTK::State stemp=s;

    //ImplicitSystemForwardEulerStep sys(2, stemp, activ);
    ImplicitSystemForwardEulerStep sys(1, stemp, excitation);

    //TODO:  Need come up with reasonable bounds
    //SimTK::Vector lower_bounds(2);
    SimTK::Vector lower_bounds(1);
    lower_bounds[0] = -SimTK::Infinity;
    //lower_bounds[1] = -SimTK::Infinity;

    //SimTK::Vector upper_bounds(2);
    SimTK::Vector upper_bounds(1);
    upper_bounds[0] = SimTK::Infinity;
    //upper_bounds[1] = SimTK::Infinity;

    sys.setParameterLimits(lower_bounds, upper_bounds);


    SimTK::Optimizer opt(sys, SimTK::InteriorPoint); //Create the optimizer

    // Specify settings for the optimizer
    opt.setConvergenceTolerance(0.1);
    opt.useNumericalGradient(true, 1e-5);
    opt.useNumericalJacobian(true);
    opt.setMaxIterations(100);
    opt.setLimitedMemoryHistory(500);

    opt.optimize(projFibVelNormGuess);  // Optimize


    double projFibVelNorm=projFibVelNormGuess[0];
    //ydot[1]=yDotInitialGuess[1];

    return projFibVelNorm;
};


SimTK::Vec3 VandenBogert2011Muscle::penMdl_calcCosPennationAngle(double projFibLenNorm, double fiberLengthNorm, bool returnJacobians) const {

    double pennAtOptFiberLength = getPennAtOptFiberLength();
    double dfiberLength_dprojFibLen = SimTK::NaN;
    double dcosPenn_dprojFibLen = SimTK::NaN;

    double cosPenn= 1.0;
    if (pennAtOptFiberLength>0.01) cosPenn = projFibLenNorm/fiberLengthNorm;

    if (returnJacobians) {
        if (pennAtOptFiberLength < 0.01) {
            // If pennation is zero, we can't do this because volume is zero,
            //      and fiberLength ~ projFiberLength
            dfiberLength_dprojFibLen = 1;
            dcosPenn_dprojFibLen = 0;
        } else {
            double b = sin(pennAtOptFiberLength);
            dfiberLength_dprojFibLen = cosPenn;
            dcosPenn_dprojFibLen = pow(b, 2) / pow(fiberLengthNorm, 3);
        }
    }

    SimTK::Vec3 cosPennAndDeriv;

    cosPennAndDeriv[0] = cosPenn;
    cosPennAndDeriv[1] = dfiberLength_dprojFibLen;
    cosPennAndDeriv[2] = dcosPenn_dprojFibLen;

    return cosPennAndDeriv;
};


double VandenBogert2011Muscle::penMdl_calcTendonLength(double projFibLengthNorm,
                                                       double muscleLength) const {
    //Returns tendonLength in meters

return muscleLength - projFibLengthNorm * getOptimalFiberLength();
};



double VandenBogert2011Muscle::penMdl_calcPennationAngularVelocity(
        double cosPennationAngle, double sinPennationAngle,
        double fiberLength, double projFiberVelocity) const
{
    double dphi = 0;
    double optimalPennationAngle = get_pennation_angle_at_optimal();

    // TODO: Note that  fiberLength*sinPennationAngle = muscle width
    // So if that is added as a fixed constant replace
    // fiberLength*sinPennationAngle and the error checks below can be moved
    // This is also the reason the pennationAngleAtOptimal is used:
    //  muscleWidthNorm = sin(pennAtOptFiberLength);

    if(get_pennation_angle_at_optimal() >= 0.01) {
        SimTK_ERRCHK_ALWAYS(fiberLength > 0,
                            "VandenBogert2011Muscle::penMdl_calcPennationAngularVelocity",
                            "Fiber length cannot be zero.");
        dphi = (fiberLength*cosPennationAngle-projFiberVelocity) /
                (fiberLength*sinPennationAngle);
    }
    return dphi;
}



double VandenBogert2011Muscle::calcFiberPotentialEnergy(double projFibLengthNorm) const {

    double optimalFiberLength=getOptimalFiberLength();
    double fiberLength=projFibLenToFiberLength(projFibLengthNorm)*optimalFiberLength;
    double fiberSlackLength = getNormFiberSlackLength()*optimalFiberLength;


    double elongationFiber = fiberLength-fiberSlackLength;


    double parrallelElementPE = -pow(elongationFiber,2) / 2;

    if (elongationFiber>0) {
        double kPEE2Norm =  getMaxIsometricForce() / (pow(getFlWidth(), 2) * optimalFiberLength);
        parrallelElementPE =
                -parrallelElementPE + (kPEE2Norm * pow(elongationFiber, 3) / 3);

    }
    return parrallelElementPE;

}


double VandenBogert2011Muscle::calcTendonPotentialEnergy(double projFibLengthNorm, double muscleLength) const {

    double projFibLength = projFibLengthNorm * getOptimalFiberLength();

    double tendonSlackLength = getTendonSlackLength();
    double elongationTendon = muscleLength - projFibLength;


    double tendonPE = -pow(elongationTendon ,2) / 2;

    if (elongationTendon>0) {

        double kSEE2Norm =  getMaxIsometricForce() / pow(getFMaxTendonStrain()*tendonSlackLength,2);

        tendonPE =
                -tendonPE + (kSEE2Norm * pow(elongationTendon, 3) / 3);

    }
    return tendonPE;

}


SimTK::Vec8 VandenBogert2011Muscle::calcFiberActiveForceVelocityMultiplier(double fiberLengtheningVelocityNorm, double cosPenn, double fiberLengthNorm, bool returnJacobians) const {

    double F2;
    double dF2_dfiberLengtheningVelocityNorm = SimTK::NaN;
    double dF2_dactiv = SimTK::NaN;
    double dF2_dprojFibVelNorm = SimTK::NaN;
    double dF2_dprojFibLenNorm = SimTK::NaN;
    double dcosPenn_dprojFibLen = SimTK::NaN;
    double dfiberLengtheningVelocityNorm_dprojFibVelNorm = SimTK::NaN;
    double dfiberLengtheningVelocityNorm_dprojFibLenNorm = SimTK::NaN;

    //AHill (dimensionless) Hill parameter of the force-velocity relationship
    double fv_AHill = getFvAHill();   //Hill parameter of the f-v relationship

    //FV_{max} (dimensionless) maximal eccentric force
    double fv_maxMultiplier = getFvmaxMultiplier();
    //Maximal eccentric force multiplier

    // Chow/Darling Vel-Activation Relationship //TODO:  Add full reference
    //double lambda = 0.5025 + 0.5341*activ;
    //double  dlambda_da = 0.5341;
    double lambda = 1;   //Turn it off for now as it seems to cause an issue
    // with negative force large concentric vel
    double dlambda_da = 0;

    if (fiberLengtheningVelocityNorm < 0) {
    //Hill's equation for concentric contraction
    // F2 = (V_{max} + V_{fiber}) / (V_{max} - V_{fiber}/a_{Hill})

        double hillDenom = (lambda * getMaxContractionVelocity() -
                            fiberLengtheningVelocityNorm / fv_AHill);
        F2 = (lambda * getMaxContractionVelocity() +
              fiberLengtheningVelocityNorm) / hillDenom;

        if (returnJacobians) {
            dF2_dfiberLengtheningVelocityNorm =
                    (1.0 + F2 / fv_AHill) / hillDenom;
            dF2_dactiv = -dlambda_da * getMaxContractionVelocity() *
                         fiberLengtheningVelocityNorm *
                         (1.0 + 1.0 / fv_AHill) / pow(hillDenom, 2);
        }
    } else {
        //Katz Curve for eccentric contraction
        // c is Katz Constant
        double c3 = getMaxContractionVelocity() * fv_AHill *
                    (fv_maxMultiplier - 1.0) /
                    (fv_AHill + 1.0); // parameter in the eccentric f-v equation
        double c = lambda * c3;
        //F2 = (g_{max} * V_{fiber} + c) / (V_{fiber} + c)
        double katzDenom = (fiberLengtheningVelocityNorm + c);
        F2 = (fv_maxMultiplier * fiberLengtheningVelocityNorm + c) / katzDenom;
        if (returnJacobians) {
            dF2_dfiberLengtheningVelocityNorm =
                    (fv_maxMultiplier - F2) / katzDenom;
            dF2_dactiv = dlambda_da * c3 * fiberLengtheningVelocityNorm *
                         (1 - fv_maxMultiplier) / pow(katzDenom, 2);
        }
    }
    if (returnJacobians) {

        double dcosPenn_dprojFibLen = pow(sin(getPennAtOptFiberLength()), 2) / pow(fiberLengthNorm, 3);

        double dfiberLengtheningVelocityNorm_dprojFibVelNorm = cosPenn;
        double dfiberLengtheningVelocityNorm_dprojFibLenNorm =
                fiberLengtheningVelocityNorm * cosPenn * dcosPenn_dprojFibLen;

        dF2_dprojFibVelNorm =
                dF2_dfiberLengtheningVelocityNorm *
                dfiberLengtheningVelocityNorm_dprojFibVelNorm;
        dF2_dprojFibLenNorm =
                dF2_dfiberLengtheningVelocityNorm *
                dfiberLengtheningVelocityNorm_dprojFibLenNorm;
    }

    SimTK::Vec8 fiberActiveForceVelocityMultiplierAndDeriv;

    fiberActiveForceVelocityMultiplierAndDeriv[0] = F2;
    fiberActiveForceVelocityMultiplierAndDeriv[1] = dF2_dfiberLengtheningVelocityNorm;
    fiberActiveForceVelocityMultiplierAndDeriv[2] = dF2_dactiv;
    fiberActiveForceVelocityMultiplierAndDeriv[3] = dF2_dprojFibVelNorm;
    fiberActiveForceVelocityMultiplierAndDeriv[4] = dF2_dprojFibLenNorm;
    fiberActiveForceVelocityMultiplierAndDeriv[5] = dcosPenn_dprojFibLen;
    fiberActiveForceVelocityMultiplierAndDeriv[6] = dfiberLengtheningVelocityNorm_dprojFibVelNorm;
    fiberActiveForceVelocityMultiplierAndDeriv[7] = dfiberLengtheningVelocityNorm_dprojFibLenNorm;

    return fiberActiveForceVelocityMultiplierAndDeriv;
}



SimTK::Vec2 VandenBogert2011Muscle::calcFiberActiveForceLengthMultiplier(double fiberLengthNorm, bool returnJacobians) const {

    double fl_width = getFlWidth();

    //---F1 is the normalized isometric force-length relationship at maximum
    //                                               activation--------------//
    double fiberExp = (fiberLengthNorm - 1.0) / fl_width;   // [dimensionless]
    double F1 = exp(-pow(fiberExp, 2));        // Gaussian force-length curve

    double dF1_dfiberLengthNorm = SimTK::NaN;

    if (returnJacobians) dF1_dfiberLengthNorm = -2.0 * fiberExp * F1 / fl_width;

    SimTK::Vec2 fiberActiveForceLengthMultiplierAndDeriv;
    fiberActiveForceLengthMultiplierAndDeriv[0] = F1;
    fiberActiveForceLengthMultiplierAndDeriv[1] = dF1_dfiberLengthNorm;

    return fiberActiveForceLengthMultiplierAndDeriv;
}


SimTK::Vec2 VandenBogert2011Muscle::calcFiberPassiveForceLengthMultiplier(double fiberLengthNorm, bool returnJacobians) const {

    //------F3 is the dimensionless fiber (PEE) force (in units of Fmax)------//
    double dF3_dfiberLengthNorm = SimTK::NaN;

    // stiffness of the linear term is 1 N/m, convert to Fmax/Lceopt units
    double kPEENorm = 1.0 / getMaxIsometricForce() * getOptimalFiberLength();;
    // elongation of fiber (PEE), relative to Lceopt
    double elongationFiberNorm = (fiberLengthNorm - getNormFiberSlackLength());
    double F3 = kPEENorm * elongationFiberNorm;
    // low stiffness linear term

    if (returnJacobians) dF3_dfiberLengthNorm = kPEENorm;

    double kPEE2Norm = 0;
    if (elongationFiberNorm > 0) {
        kPEE2Norm = 1 / pow(getFlWidth(), 2);  //Fiber (PEE) quadratic stiffness,
        //          so Fpee = Fmax when Lce = Lce*(1+W)
        //add quadratic term for positive elongation
        F3 = F3 + kPEE2Norm * pow(elongationFiberNorm, 2);
        if (returnJacobians) {
            dF3_dfiberLengthNorm =
                    dF3_dfiberLengthNorm + 2 * kPEE2Norm * elongationFiberNorm;

        }
    }

    SimTK::Vec2 fiberPassiveForceLengthMultiplierAndDeriv;
    fiberPassiveForceLengthMultiplierAndDeriv[0]=F3;
    fiberPassiveForceLengthMultiplierAndDeriv[1]=dF3_dfiberLengthNorm;

    return fiberPassiveForceLengthMultiplierAndDeriv;
}


SimTK::Vec4 VandenBogert2011Muscle::calcTendonForceMultiplier(double muscleLength, double projFibLenNorm, bool returnJacobians) const {

// elongationTendon in [m]
//--------F4 is the dimensionless SEE force (in units of Fmax)----------//

//TODO: kSEE and kSEE2 could probably become derived parameters

// stiffness of the linear term is 1 N/m, convert to Fmax/m (so normalized
// by Fmax)
    double kSEE = 1.0 / getMaxIsometricForce();
    double optFiberLength = getOptimalFiberLength();
    double elongationTendon = muscleLength - projFibLenNorm * optFiberLength; //m

//  low stiffness linear term
    double F4 = kSEE * elongationTendon;
    double dF4_dmuscleLength = SimTK::NaN;
    double dF4_dprojFibLenNorm = SimTK::NaN;
    double dF4_dtendonLength = SimTK::NaN;
    if (returnJacobians) {
        dF4_dprojFibLenNorm = -kSEE * optFiberLength;
        dF4_dmuscleLength = kSEE;
        dF4_dtendonLength = kSEE;  //In N/m units
    }

    double kSEE2 = SimTK::NaN;
    if (elongationTendon > 0) {
// Tendon (SEE) quadratic stiffness, so Fsee = Fmax at strain of umax
// This is normalized by Fmax
        kSEE2 = 1 / pow(getTendonSlackLength() * getFMaxTendonStrain(), 2);


// add quadratic term for positive deformation
        F4 = F4 + (kSEE2 * pow(elongationTendon, 2));
        if (returnJacobians) {
            dF4_dprojFibLenNorm = dF4_dprojFibLenNorm - 2 * kSEE2 *
                                                        optFiberLength *
                                                        elongationTendon;
            dF4_dmuscleLength = dF4_dmuscleLength + 2 * kSEE2 *
                                                    elongationTendon;
            dF4_dtendonLength = dF4_dtendonLength  + 2*kSEE2*elongationTendon;
        }
    }


    SimTK::Vec4 tendonForceMultiplierAndDeriv;
    tendonForceMultiplierAndDeriv[0]=F4;
    tendonForceMultiplierAndDeriv[1]=dF4_dprojFibLenNorm;
    tendonForceMultiplierAndDeriv[2]=dF4_dmuscleLength;
    tendonForceMultiplierAndDeriv[3]=dF4_dtendonLength;

    return tendonForceMultiplierAndDeriv;
}


SimTK::Vec2 VandenBogert2011Muscle::calcFiberDampingMultiplier(double projFibVelNorm, bool returnJacobians) const {

    double dampingCoeff = getDampingCoefficient();
    double F5 = dampingCoeff * projFibVelNorm ;
    double dF5_dprojFibVelNorm  = dampingCoeff;


    SimTK::Vec2 fiberDampingMultiplierAndDeriv;
    fiberDampingMultiplierAndDeriv[0] = F5;
    fiberDampingMultiplierAndDeriv[1] = dF5_dprojFibVelNorm;

    return fiberDampingMultiplierAndDeriv;
}
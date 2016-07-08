/* -------------------------------------------------------------------------- *
 *                      OpenSim:  VandenBogert.cpp                            *
 * -------------------------------------------------------------------------- */


//=============================================================================
// INCLUDES
//=============================================================================
#include "VandenBogert2011Muscle.h"

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
}

//_____________________________________________________________________________
/*
 * Constructor.
 */
VandenBogert2011Muscle::VandenBogert2011Muscle(const std::string &name)
{
    constructProperties();
}

//_____________________________________________________________________________
/*
 * Construct and initialize properties.
 * All properties are added to the property set. Once added, they can be
 * read in and written to files.
 */
void VandenBogert2011Muscle::constructProperties()
{
    setAuthors("A. van den Bogert, B. Humphreys");
    //constructProperty_max_isometric_force(1000);
    constructProperty_fMaxTendonStrain(0.04);
    constructProperty_fl_width(0.63);
    constructProperty_fv_AHill(0.25);
    constructProperty_fv_maxMultiplier(1.4);
    //constructProperty_optimal_fiber_length(0.1);
    constructProperty_dampingCoefficient(0.01);
    constructProperty_normFiberSlackLength(1.0);
    //constructProperty_tendon_slack_length(0.2);
    constructProperty_activTimeConstant(0.01);
    constructProperty_deactivTimeConstant(0.04);
    constructProperty_pennAtOptFiberLength(0.0);
}

// Define new states and their derivatives in the underlying system
void VandenBogert2011Muscle::extendAddToSystem(SimTK::MultibodySystem& system) const
{

    // No States to add, yet
    /* Allow Millard2012EquilibriumMuscle to add its states, before extending
    Super::extendAddToSystem(system);

    // Now add the states necessary to implement the fatigable behavior
    addStateVariable("target_activation");
    addStateVariable("active_motor_units");
    addStateVariable("fatigued_motor_units");
    // and their corresponding derivatives
    addCacheVariable("target_activation_deriv", 0.0, SimTK::Stage::Dynamics);
    addCacheVariable("active_motor_units_deriv", 0.0, SimTK::Stage::Dynamics);
    addCacheVariable("fatigued_motor_units_deriv", 0.0, SimTK::Stage::Dynamics);*/
}

void VandenBogert2011Muscle::extendInitStateFromProperties(SimTK::State& s) const
{
    // No States Yet
    /*Super::extendInitStateFromProperties(s);
    setTargetActivation(s, getDefaultTargetActivation());
    setActiveMotorUnits(s, getDefaultActiveMotorUnits());
    setFatiguedMotorUnits(s, getDefaultFatiguedMotorUnits());*/
}

void VandenBogert2011Muscle::extendSetPropertiesFromState(const SimTK::State& s)
{
    // No States Yet
    /*Super::extendSetPropertiesFromState(s);
    setDefaultTargetActivation(getTargetActivation(s));
    setDefaultActiveMotorUnits(getActiveMotorUnits(s));
    setDefaultFatiguedMotorUnits(getFatiguedMotorUnits(s));*/
}

//--------------------------------------------------------------------------
// GET & SET Properties
//--------------------------------------------------------------------------
void VandenBogert2011Muscle::setfMaxTendonStrain(double fMaxTendonStrain) {
    setfMaxTendonStrain(fMaxTendonStrain);
    }
double VandenBogert2011Muscle::getfMaxTendonStrain() const
    { return getfMaxTendonStrain(); }

void VandenBogert2011Muscle::setfl_width(double fl_width)
    { setfl_width(fl_width); }
double VandenBogert2011Muscle::getfl_width() const
    { return getfl_width(); }

void VandenBogert2011Muscle::setfv_AHill(double fv_AHill)
    { setfv_AHill(fv_AHill); }
double VandenBogert2011Muscle::getfv_AHill() const
    { return getfv_AHill(); }

void VandenBogert2011Muscle::setfv_maxMultiplier(double fv_maxMultiplier)
    { setfv_maxMultiplier(fv_maxMultiplier); }
double VandenBogert2011Muscle::getfv_maxMultiplier() const
    { return getfv_maxMultiplier(); }

void VandenBogert2011Muscle::setdampingCoefficient(double dampingCoefficient)
    { setdampingCoefficient(dampingCoefficient); }
double VandenBogert2011Muscle::getdampingCoefficient() const
    { return getdampingCoefficient(); }

void VandenBogert2011Muscle::setnormFiberSlackLength(double normFiberSlackLength)
    { setnormFiberSlackLength(normFiberSlackLength); }
double VandenBogert2011Muscle::getnormFiberSlackLength() const
    { return getnormFiberSlackLength(); }

void VandenBogert2011Muscle::setactivTimeConstant(double activTimeConstant)
    { setactivTimeConstant(activTimeConstant); }
double VandenBogert2011Muscle::getactivTimeConstant() const
    { return getactivTimeConstant(); }

void VandenBogert2011Muscle::setdeactivTimeConstant(double deactivTimeConstant)
    { setdeactivTimeConstant(deactivTimeConstant); }
double VandenBogert2011Muscle::getdeactivTimeConstant() const
    { return getdeactivTimeConstant(); }

void VandenBogert2011Muscle::setpennAtOptFiberLength(double pennAtOptFiberLength)
{ setpennAtOptFiberLength(pennAtOptFiberLength); }
double VandenBogert2011Muscle::getpennAtOptFiberLength() const
{ return getpennAtOptFiberLength(); }


//--------------------------------------------------------------------------
// GET & SET States and their derivatives
//--------------------------------------------------------------------------

//No states yet
/*double FatigableMuscle::getTargetActivation(const SimTK::State& s) const
{   return getStateVariableValue(s, "target_activation"); }

void FatigableMuscle::setTargetActivation(SimTK::State& s,
                                          double fatiguedAct) const
{   setStateVariableValue(s, "target_activation", fatiguedAct); }

double FatigableMuscle::getTargetActivationDeriv(const SimTK::State& s) const
{   return getStateVariableDerivativeValue(s, "target_activation"); }

void FatigableMuscle::setTargetActivationDeriv(const SimTK::State& s,
                                               double fatiguedActDeriv) const
{   setStateVariableDerivativeValue(s, "target_activation", fatiguedActDeriv); }

double FatigableMuscle::getActiveMotorUnits(const SimTK::State& s) const
{   return getStateVariableValue(s, "active_motor_units"); }

void FatigableMuscle::setActiveMotorUnits(SimTK::State& s,
                                          double activeMotorUnits) const
{   setStateVariableValue(s, "active_motor_units", activeMotorUnits); }

double FatigableMuscle::getActiveMotorUnitsDeriv(const SimTK::State& s) const
{   return getStateVariableDerivativeValue(s, "active_motor_units"); }

void FatigableMuscle::setActiveMotorUnitsDeriv(const SimTK::State& s,
                                               double activeMotorUnitsDeriv) const
{   setStateVariableDerivativeValue(s, "active_motor_units", activeMotorUnitsDeriv); }

double FatigableMuscle::getFatiguedMotorUnits(const SimTK::State& s) const
{   return getStateVariableValue(s, "fatigued_motor_units"); }

void FatigableMuscle::setFatiguedMotorUnits(SimTK::State& s,
                                            double fatiguedMotorUnits) const
{   setStateVariableValue(s, "fatigued_motor_units", fatiguedMotorUnits); }

double FatigableMuscle::getFatiguedMotorUnitsDeriv(const SimTK::State& s) const
{    return getStateVariableDerivativeValue(s, "fatigued_motor_units"); }

void FatigableMuscle::setFatiguedMotorUnitsDeriv(const SimTK::State& s,
                                                 double fatiguedMotorUnitsDeriv) const
{   setStateVariableDerivativeValue(s, "fatigued_motor_units", fatiguedMotorUnitsDeriv);
 }*/





//==============================================================================
// Muscle.h Interface
//==============================================================================


double  VandenBogert2011Muscle::computeActuation(const SimTK::State& s) const
{return( 0.0);}

void VandenBogert2011Muscle::setActivation(SimTK::State& s,double activation) const
{}

void VandenBogert2011Muscle::computeInitialFiberEquilibrium(SimTK::State& s) const
{}




//=============================================================================
// COMPUTATION
//=============================================================================
//_____________________________________________________________________________
/**
 * Compute the derivatives of the muscle states.
 *
 * @param s  system state
 */


//array<double,3> calcImplicitResidual(const SimTK::State& s) const

//SimTK::Vec3 VandenBogert2011Muscle::calcImplicitResidual(double Lm, double Lce, double a, double Lcedot, double adot, double u, int returnJacobians) const {


VandenBogert2011Muscle::ImplicitResults VandenBogert2011Muscle::calcImplicitResidual(double muslceLength, double projFibLen, double activ, double projFibVel, double activ_dot, double u, int returnJacobians) const {

    // Later when using state variables:
    //State Variables
    //double Lm = getLength(s);
    //double Lce = getNormalizedFiberLength(s);  //Ton: Lce is dimensionless, it is the muscle fiber length divided by Lceopt
    //double a = getActivation(s);
    //double Lcedot = getNormalizedFiberVelocity(s);

    //TODO: Match symbols to doxygen diagram
    // -------------------------Parameters----------------------------//
    //F_{iso} Maximum isometric force that the fibers can generate
    double maxIsoForce = getMaxIsometricForce();

    //u_{max} (dimensionless) strain in the series elastic element at load of maxIsometricForce
    double fMaxTendonStrain = getfMaxTendonStrain(); //Strain in the series elastic element at load of Fmax

    //W (dimensionless) width parameter of the force-length relationship of the muscle fiber
    double fl_width = getfl_width();

    //AHill (dimensionless) Hill parameter of the force-velocity relationship
    double fv_AHill= getfv_AHill();        //Hill parameter of the force-velocity relationship

    //FV_{max} (dimensionless) maximal eccentric force
    double fv_maxMultiplier = getfv_maxMultiplier(); //Maximal eccentric force multiplier

    // L_{opt}  (m) Optimal Length of Contractile Element;
    double optFiberLength = getOptimalFiberLength();  //Optimal Length of Contractile Element

    //b (s/m) damping coefficient of damper parallel to the fiber (normalized to maxIsometricForce)
    double dampingCoeff = getdampingCoefficient();

    //L_{slack,fiber} Slack length of the parallel elastic element, divided by Lceopt
    double fiberSlackLength= getnormFiberSlackLength();

    //L_{slack,tendon} (dimensionless) slack length of the tendon
    double tendonSlackLength = get_tendon_slack_length();

    //T_{act} (s) Activation time
    double activTimeConstant = getactivTimeConstant();

    //T_{deact} (s) Deactivation time
    double deactivTimeConstant = getdeactivTimeConstant();

    //phi_{opt} pennation at Optimal Fiber Length
    double pennAtOptFiberLength = getpennAtOptFiberLength();


    // constants derived from the muscle parameters
    double Vmax = 10 * optFiberLength;                 //Maximum shortening velocity (m/s) is 10 fiber lengths per sec


    double kPEE2 = 1 / pow(fl_width, 2);                       // Fiber (PEE) quadratic stiffness, so Fpee = Fmax when Lce = Lce*(1+W)
    double kSEE2 = 1 / pow(fiberSlackLength * fMaxTendonStrain, 2);  // Tendon (SEE) quadratic stiffness, so Fsee = Fmax at strain of umax


    // Jacobian Matrices
    SimTK::Mat23 df_dy;
    SimTK::Mat23 df_dydot;


    //-------------------Convert projFiberLength & Velocity to fiberLength & Velocity------/
    //TODO: May want to make this into a seperate function
    //TODO:  Add symbolic equations comments
    double cosPenn;
    double dcosPenn_dprojFibLen;
    double fiberLength;
    double dfiberLength_dprojFibLen;
    if (pennAtOptFiberLength<0.01)  {
        // If pennation is zero, we can't do this because volume is zero, and fiberLength ~ projFiberLength
        cosPenn=1.0;
        fiberLength=projFibLen;
        dfiberLength_dprojFibLen=1;
        dcosPenn_dprojFibLen=0;
    }
    else {
        double b=sin(pennAtOptFiberLength);
        fiberLength=sqrt(pow(projFibLen,2) + pow(b,2));
        cosPenn = projFibLen/fiberLength;
        dfiberLength_dprojFibLen=cosPenn;
        dcosPenn_dprojFibLen=pow(b,2) / pow(fiberLength,3);

    }

    // Compute fiberVelocity and its derivatives wrt projFibLen and projFibVel
    double fiberVelocity = projFibVel*cosPenn;
    double dfiberVelocity_dprojFibVel = cosPenn;
    double dfiberVelocity_dprojFibLen = projFibVel * dcosPenn_dprojFibLen;



    //----------------F1 is the normalized isometric force-length relationship at maximum activation--------------//
    double fiberExp = (fiberLength - 1.0) / fl_width;   // [dimensionless]
    double F1 = exp(pow(-fiberExp, 2));        // Gaussian force-length curve

    double dF1_dfiberLength = 0;
    if (returnJacobians) {
        double dF1_dfiberLength = -2.0 * fiberExp * F1 / fl_width;
        double dF1_dprojFibLen = dF1_dfiberLength * dfiberLength_dprojFibLen;
    }



    //----------------- F2 is the dimensionless force-velocity relationship --------------------------------//
    double F2;
    double dF2_dfiberVelocity;
    double dF2_dactiv;
    double dF2_dprojFibVel;
    double dF2_dprojFibLen;

    // Chow/Darling Vel-Activation Relationship //TODO:  Add full reference
    //double lambda = 0.5025 + 0.5341*activ;
    //double  dlambda_da = 0.5341;
    double lambda = 1;   //Turn it off for now as it seems to cause an issue with negative force large concentric vel
    double dlambda_da =0;

    if (fiberVelocity < 0) {
        //Hill's equation for concentric contraction
        // F2 = (V_{max} + V_{fiber}) / (V_{max} - V_{fiber}/a_{Hill})
        double hillDenom = (lambda*Vmax - fiberVelocity/fv_AHill);
        F2 = (lambda*Vmax + fiberVelocity) / hillDenom;
        if (returnJacobians) {
            dF2_dfiberVelocity  = (1.0 + F2 / fv_AHill) / hillDenom;
            dF2_dactiv = - dlambda_da * Vmax * fiberVelocity * (1 + 1/fv_AHill) / pow(hillDenom,2);
        }
    }
    else {
        //Katz Curve for eccentric contraction
        // c is Katz Constant
        double c3 = Vmax * fv_AHill * (fv_maxMultiplier - 1.0) / (fv_AHill + 1.0); // parameter in the eccentric force-velocity equation
        double c = lambda*c3;
        //F2 = (g_{max} * V_{fiber} + c) / (V_{fiber} + c)
        double katzDenom = (fiberVelocity  + c);
        F2 = (fv_maxMultiplier * fiberVelocity  + c) / katzDenom ;
        if (returnJacobians) {
            dF2_dfiberVelocity = (fv_maxMultiplier - F2) / katzDenom ;
            dF2_dactiv = dlambda_da * c3 * fiberVelocity * (1-fv_maxMultiplier) / pow(katzDenom,2);
        }
    }
    if (returnJacobians){
        dF2_dprojFibVel = dF2_dfiberVelocity * dfiberVelocity_dprojFibVel;
        dF2_dprojFibLen = dF2_dfiberVelocity * dfiberVelocity_dprojFibLen;
    }


    //----------F3 is the dimensionless fiber (PEE) force (in units of Fmax)-------------//
    double dF3_dprojFibLen;
    double kPEE = 1.0 / maxIsoForce * optFiberLength;      // stiffness of the linear term is 1 N/m, convert to Fmax/Lceopt units
    double elongationFiber = (fiberLength - fiberSlackLength);  // elongation of fiber (PEE), relative to Lceopt
    double F3 = kPEE * elongationFiber;                // low stiffness linear term
    double dF3_dfiberLength = kPEE;
    if (elongationFiber > 0) {
        //add quadratic term for positive elongation
        F3 = F3 + kPEE2 * pow(elongationFiber, 2);
        if (returnJacobians) {
            dF3_dfiberLength = dF3_dfiberLength + 2 * kPEE2 * elongationFiber;
        }
    }
    if (returnJacobians) {
        dF3_dprojFibLen = dF3_dfiberLength * dfiberLength_dprojFibLen;
    }


    //---------------F4 is the dimensionless SEE force (in units of Fmax)----------//
    double kSEE = 1.0 / maxIsoForce;                  // stiffness of the linear term is 1 N/m, convert to Fmax/m
    double elongationTendon = muslceLength - fiberLength * optFiberLength - tendonSlackLength;  // elongation of tendon (SEE), in meters
    double F4 = kSEE * elongationTendon;              //  low stiffness linear term
    double dF4_dfiberLength;
    double dF4_dmuscleLength;
    double dF4_dprojFibLen;
    if (returnJacobians) {
         double dF4_dprojFibLen = -kSEE * optFiberLength;
         dF4_dmuscleLength = kSEE;
    }
    if (elongationTendon > 0) {
        // add quadratic term for positive deformation
        F4 = F4 + kSEE2 * pow(elongationTendon, 2);
        if (returnJacobians) {
             dF4_dprojFibLen = dF4_dprojFibLen - 2 * kSEE2 * optFiberLength* elongationTendon;
             dF4_dmuscleLength = dF4_dmuscleLength + 2 * kSEE2 * elongationTendon;
        }
    }



    //-------------- F5 is viscous damping parallel to the CE (0.001 of Fmax at 1 Lceopt/s) to  -----------//
    // ensure that df/dLcedot is never zero
    double F5 = dampingCoeff * projFibVel ;   //TODO:  Should this be: dampingCoeff * fiberVelocity  (not proj)?
    double dF5_dprojFibVel  = dampingCoeff;



    // ----------------------Calculate the Muscle Force Residual ------------------------- //
    //The muscle dynamics equation: f = Fsee - (a*Fce - Fpee)*cos(Penn) - Fdamping = 0
    double fRes = F4 - (activ * F1 * F2 + F3)*cosPenn - F5;



    //---------------------Assemble Force Jacobian---------------------------//
    if (returnJacobians) {

        double dfRes_dactiv = -(F1*F2 + activ*F1*dF2_dactiv )*cosPenn;
        double dfRes_dprojFibLength = dF4_dprojFibLen - (activ*(dF1_dfiberLength*F2 + F1*dF2_dprojFibLen) + dF3_dprojFibLen) * cosPenn - (activ*F1*F2 + F3)*dcosPenn_dprojFibLen;
        double dfRes_dprojFibVel = - activ*F1*dF2_dprojFibVel - dF5_dprojFibVel;
        double df_dmuscleLength = dF4_dmuscleLength;

        df_dy[0][0] = df_dmuscleLength;
        df_dy[0][1] = dfRes_dprojFibLength;
        df_dy[0][1] = dfRes_dactiv;

        df_dy[1][0] = 0;
        df_dy[1][1] = 0;
        df_dy[1][2] = dfRes_dprojFibVel;

    }



    // ---------------- Force in tendon (SEE) is maxIsoForce*F4 ---------------------------- //
    double Fsee = maxIsoForce * F4;

    // TODO: I don't think we need to implement dFsee_.....
    /*if (nargout > 1)
        dFsee_dLm  = Fmax * dF4_dLm;
    dFsee_dLce = Fmax * dF4_dLce;
    dFsee_dy   = [0;0;(-d/L)*dFsee_dLm; 0;dFsee_dLce;0;0];  % Fsee is a function of y(3) & y(5)
    end*/



    //----------------------Activation dynamics equation-------------------//
    double df_du = 0;
    double activationResidual = activ_dot - (u - activ) * (u / activTimeConstant + (1 - u) / deactivTimeConstant );
    if (returnJacobians) {
        df_dy[1][2] = (u / activTimeConstant + (1 - u) / deactivTimeConstant );
        double df_du = -(u / activTimeConstant + (1 - u) / deactivTimeConstant ) - (u - activ) * (1 / activTimeConstant - 1 / deactivTimeConstant );
        df_dydot[1][2] = 1;
    }



    VandenBogert2011Muscle::ImplicitResults Results;

    Results.forceResidual=fRes;
    Results.activResidual=activationResidual;
    Results.forceTendon=Fsee;
    Results.df_dy=df_dy;
    Results.df_dydot=df_dydot;
    Results.df_du=df_du;
    Results.F1=F1;    //Output these for troubleshooting
    Results.F2=F2;
    Results.F3=F3;
    Results.F4=F4;
    Results.F5=F5;



return Results;

}

/* Author: Matthew Millard
/*
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include <OpenSim/Actuators/ActiveForceLengthCurve.h>
#include <OpenSim/Actuators/ForceVelocityCurve.h>
#include <OpenSim/Actuators/ForceVelocityInverseCurve.h>
#include <OpenSim/Actuators/TendonForceLengthCurve.h>
#include <OpenSim/Actuators/FiberForceLengthCurve.h>
#include <OpenSim/Actuators/FiberCompressiveForceLengthCurve.h>
#include <OpenSim/Actuators/FiberCompressiveForceCosPennationCurve.h>

#include <SimTKsimbody.h>
#include <ctime>
#include <string>
#include <stdio.h>

using namespace std;
using namespace OpenSim;
using namespace SimTK;

void testActiveForceLengthCurve();
void testForceVelocityCurve();
void testForceVelocityInverseCurve();
void testTendonForceLengthCurve();
void testFiberForceLengthCurve();
void testFiberCompressiveForceLengthCurve();
void testFiberCompressiveForceCosPennationCurve();

int main(int argc, char* argv[])
{
	

	try {
            SimTK_START_TEST("Testing Serializable Curves");
            testActiveForceLengthCurve();
            testForceVelocityCurve();
            testForceVelocityInverseCurve();
            testTendonForceLengthCurve();
            testFiberForceLengthCurve();
            testFiberCompressiveForceLengthCurve();
            testFiberCompressiveForceCosPennationCurve();
            SimTK_END_TEST();
    }
    catch (OpenSim::Exception ex)
    {
        cout << ex.getMessage() << endl;
		cin.get();
        return 1;
    }
    catch (const std::exception& ex)
    {
        cout << ex.what() << endl;
		cin.get();		
		return 1;
    }
    catch (...)
    {
        cout << "UNRECOGNIZED EXCEPTION" << endl;
		cin.get();
        return 1;
    }

	

    cout << "\n Serializable Curve Testing completed successfully.\n";
	return 0;
}

void testActiveForceLengthCurve()
{

        //cout <<"**************************************************"<<endl;
        cout <<"________________________________________________________"<<endl;
        cout <<"1. Testing: ActiveForceLengthCurve "<<endl;       
        cout <<"________________________________________________________"<<endl;

        cout <<"    a. default construction" <<endl;
        ActiveForceLengthCurve falCurve1;
        //falCurve1.setName("default_ActiveForceLengthCurve");
        falCurve1.print("default_ActiveForceLengthCurve.xml");

        cout <<"    b. serialization & deserialization" <<endl;
        ActiveForceLengthCurve falCurve2;
        falCurve2.setMaxActiveFiberLength(2);
        falCurve2.setTransitionFiberLength(0.8);
        falCurve2.setMinActiveFiberLength(0);
        falCurve2.setMinValue(0.3);
        falCurve2.setShallowAscendingSlope(0.5);


        //These next few lines are just to read the object in, and repopulate
        //falCurve2 with the properties from the file ... and its a little 
        //awkward to use.

        //cout << "b.*Uncomment, test makeObjectFromFile once in OpenSim"<<endl;        
        
        
        Object* tmpObj = Object::
                       makeObjectFromFile("default_ActiveForceLengthCurve.xml");
        falCurve2 = *dynamic_cast<ActiveForceLengthCurve*>(tmpObj);
        delete tmpObj;
        
        SimTK_TEST(falCurve2 == falCurve1);        
        remove("default_ActiveForceLengthCurve.xml");  
        

        falCurve2.setMaxActiveFiberLength(2);
        falCurve2.setTransitionFiberLength(0.8);
        falCurve2.setMinActiveFiberLength(0);
        falCurve2.setMinValue(0.3);
        falCurve2.setShallowAscendingSlope(0.5);

        cout <<"    c. assignment operator" <<endl;
        falCurve2=falCurve1;
        
        
        SimTK_TEST(falCurve1==falCurve2);

        falCurve2.setMaxActiveFiberLength(2);
        falCurve2.setTransitionFiberLength(0.8);
        falCurve2.setMinActiveFiberLength(0);
        falCurve2.setMinValue(0.3);
        falCurve2.setShallowAscendingSlope(0.5);

        cout <<"    d. copy constructor" <<endl;
        ActiveForceLengthCurve falCurve2p5(falCurve2);
        SimTK_TEST(falCurve2==falCurve2p5);

        cout << "*Passed: default construction, limited serialization" << endl;
        cout << "         assignment operator, copy constructor" << endl;

        //cout <<"**************************************************"<<endl;
        cout <<"2. Testing API constructor" << endl;
        ActiveForceLengthCurve falCurve3(0.5, 0.75,1.5,0.75,0.01,"testMuscle");
        double falVal  = falCurve3.calcValue(1.0);
        double dfalVal = falCurve3.calcDerivative(1.0,1);
        cout << "Passed: Testing API constructor" << endl;

        //cout <<"**************************************************"<<endl;
        cout <<"3. Testing get/set methods:" << endl;

        falCurve2.setMinActiveFiberLength(0);
        falCurve2.setTransitionFiberLength(0.8);
        falCurve2.setMaxActiveFiberLength(2);        
        falCurve2.setMinValue(0.3);
        falCurve2.setShallowAscendingSlope(0.5);

        SimTK_TEST(falCurve2.getMinActiveFiberLength() == 0.0);
        SimTK_TEST(falCurve2.getTransitionFiberLength()== 0.8);
        SimTK_TEST(falCurve2.getMaxActiveFiberLength() == 2.0);
        SimTK_TEST(falCurve2.getMinValue() == 0.3);
        SimTK_TEST(falCurve2.getShallowAscendingSlope() == 0.5);

        cout << "Passed: Testing get/set methods" << endl;

        //====================================================================
        double p1 = 0.4;
        double p2 = 0.75;
        double p3 = 1.6;
        double p4 = 0.75;
        double p5 = 0.05;


        printf("4. Testing default property values: \n\t%f,\n\t%f,\n\t%f,"
               "\n\t%f,\n\t%f\n",p1,p2,p3,p4,p5);      
        
        ActiveForceLengthCurve falCurve4;
        falCurve4.setName("falCurve");

        SimTK_TEST(falCurve4.getMinActiveFiberLength()  == p1);
        SimTK_TEST(falCurve4.getTransitionFiberLength() == p2);
        SimTK_TEST(falCurve4.getMaxActiveFiberLength()  == p3);
        SimTK_TEST(falCurve4.getShallowAscendingSlope() == p4);
        SimTK_TEST(falCurve4.getMinValue()              == p5);
        cout << "Passed" << endl;
        //====================================================================

        cout <<"5. Testing Services for connectivity:" << endl;        
        

        cout <<"    a. calcValue" << endl;
            double tol = sqrt(SimTK::Eps);
            double value = falCurve4.calcValue(1.0);
            SimTK_TEST_EQ_TOL(value, 1, tol);
        cout <<"    b. calcDerivative" << endl;
            double dvalue= falCurve4.calcDerivative(1.0,1);
            SimTK_TEST_EQ_TOL(dvalue, 0, tol);

        cout <<"    c. getCurveDomain" << endl;
            SimTK::Vec2 tmp = falCurve4.getCurveDomain();
            SimTK_TEST(tmp(0) == p1 &&
                       tmp(1) == p3);

        cout <<"    d. printMuscleCurveToCSVFile" << endl;
            falCurve4.printMuscleCurveToCSVFile("");
            std::string fname = falCurve4.getName();
            fname.append(".csv");
            remove(fname.c_str());

        cout << "Passed: Testing Services for connectivity" << endl;                            

        //cout <<"**************************************************"<<endl;
        cout <<"Service correctness is tested by underlying utility class"<<endl;
        cout <<"MuscleCurveFunction, and MuscleCurveFunctionFactory"<<endl;
        //cout <<"**************************************************"<<endl;

        //cout <<"**************************************************"<<endl;
        cout <<"________________________________________________________"<<endl;
        cout <<"          TESTING ActiveForceLengthCurve          "<<endl;
        cout <<"                    COMPLETED                     "<<endl;
        cout <<"________________________________________________________"<<endl;
        //cout <<"**************************************************"<<endl;

}

void testForceVelocityCurve()
{

        cout <<"________________________________________________________"<<endl;
        cout <<"1. Testing ForceVelocityCurve"<<endl;       
        cout <<"________________________________________________________"<<endl;


        cout <<"    a. default construction" <<endl;
        ForceVelocityCurve fvCurve1;
        fvCurve1.print("default_ForceVelocityCurve.xml");

        cout <<"    b. serialization & deserialization" <<endl;
        ForceVelocityCurve fvCurve2;
        //change all of the properties to something other than the default
        fvCurve2.setConcentricCurviness(0.5);
        fvCurve2.setConcentricMinSlope(0);
        fvCurve2.setEccentricCurviness(0.5);
        fvCurve2.setEccentricMinSlope(0);
        fvCurve2.setMaxEccentricVelocityForceMultiplier(2.0);
        fvCurve2.setIsometricMaxSlope(10);


        //These next few lines are just to read the object in, and repopulate
        //fvCurve2 with the properties from the file ... and its a little 
        //awkward to use.

        //cout << "b.*Uncomment, test makeObjectFromFile once in OpenSim"<<endl;        
        
        
        Object* tmpObj = Object::          
                       makeObjectFromFile("default_ForceVelocityCurve.xml");
        fvCurve2 = *dynamic_cast<ForceVelocityCurve*>(tmpObj);        
        delete tmpObj;
        SimTK_TEST(fvCurve2 == fvCurve1);       
        remove("default_ForceVelocityCurve.xml");
        

        fvCurve2.setConcentricCurviness(0.5);
        fvCurve2.setConcentricMinSlope(0);
        fvCurve2.setEccentricCurviness(0.5);
        fvCurve2.setEccentricMinSlope(0);
        fvCurve2.setMaxEccentricVelocityForceMultiplier(2.0);
        fvCurve2.setIsometricMaxSlope(10);

        cout <<"    c. assignment operator" <<endl;
        fvCurve2=fvCurve1;
        
        
        SimTK_TEST(fvCurve1==fvCurve2);

        fvCurve2.setConcentricCurviness(0.5);
        fvCurve2.setConcentricMinSlope(0);
        fvCurve2.setEccentricCurviness(0.5);
        fvCurve2.setEccentricMinSlope(0);
        fvCurve2.setMaxEccentricVelocityForceMultiplier(2.0);
        fvCurve2.setIsometricMaxSlope(10);

        cout <<"    d. copy constructor" <<endl;
        ForceVelocityCurve fvCurve2p5(fvCurve2);
        SimTK_TEST(fvCurve2==fvCurve2p5);

        cout << "*Passed: default construction, limited serialization" << endl;
        cout << "         assignment operator, copy constructor" << endl;

        //====================================================================

        cout <<"2. Testing API constructor" << endl;
        ForceVelocityCurve fvCurve3(0,5,0,1.8,0.1,0.75,"testMuscle");
        double falVal  = fvCurve3.calcValue(1.0);
        double dfalVal = fvCurve3.calcDerivative(1.0,1);
        cout << "Passed: Testing API constructor" << endl;

        //====================================================================

        cout <<"3. Testing get/set methods:" << endl;

        fvCurve2.setConcentricCurviness(0.5);
        fvCurve2.setConcentricMinSlope(0);
        fvCurve2.setEccentricCurviness(0.6);
        fvCurve2.setEccentricMinSlope(0.1);
        fvCurve2.setMaxEccentricVelocityForceMultiplier(2.0);
        fvCurve2.setIsometricMaxSlope(10);

        SimTK_TEST(fvCurve2.getConcentricCurviness()                    == 0.5);
        SimTK_TEST(fvCurve2.getConcentricMinSlope()                     == 0  );
        SimTK_TEST(fvCurve2.getEccentricCurviness()                     == 0.6);
        SimTK_TEST(fvCurve2.getEccentricMinSlope()                      == 0.1);
        SimTK_TEST(fvCurve2.getMaxEccentricVelocityForceMultiplier()    == 2.0);
        SimTK_TEST(fvCurve2.getIsometricMaxSlope()                      ==  10);

        cout << "Passed: Testing get/set methods" << endl;

        //====================================================================
        double p1 = 0.1;
        double p2 = 5;
        double p3 = 0.1;
        double p4 = 1.8;
        double p5 = 0.1;
        double p6 = 0.75;

        printf("4. Testing default property values: "
               "\n\t%f,\n\t%f,\n\t%f,\n\t%f,\n\t%f,\n\t%f\n"
            ,p1,p2,p3,p4,p5,p6);      
        ForceVelocityCurve fvCurve4;
        fvCurve4.setName("fvCurve");

        SimTK_TEST(fvCurve4.getConcentricMinSlope()     == p1);
        SimTK_TEST(fvCurve4.getIsometricMaxSlope()      == p2);
        SimTK_TEST(fvCurve4.getEccentricMinSlope()      == p3);
        SimTK_TEST(fvCurve4.getMaxEccentricVelocityForceMultiplier() == p4);
        SimTK_TEST(fvCurve4.getConcentricCurviness()    == p5);
        SimTK_TEST(fvCurve4.getEccentricCurviness()     == p6);
        cout << "Passed" << endl;
        //====================================================================

        cout <<"5. Testing Services for connectivity:" << endl;        

        cout <<"    a. calcValue" << endl;
            double tol = sqrt(SimTK::Eps);
            double value = fvCurve4.calcValue(0);
            SimTK_TEST_EQ_TOL(value,1, tol);
        cout <<"    b. calcDerivative" << endl;
            double dvalue= fvCurve4.calcDerivative(0,1);
            SimTK_TEST_EQ_TOL(dvalue, 5, tol);
            dvalue= fvCurve4.calcDerivative(0,2);
            SimTK_TEST_EQ_TOL(dvalue, 0, tol);

        cout <<"    c. getCurveDomain" << endl;
            SimTK::Vec2 tmp = fvCurve4.getCurveDomain();
            SimTK_TEST(tmp(0) == -1.0 &&
                       tmp(1) == 1.0);

        cout <<"    d. printMuscleCurveToCSVFile" << endl;
            fvCurve4.setConcentricCurviness(0.5);
            fvCurve4.setEccentricCurviness(1.0);
            fvCurve4.printMuscleCurveToCSVFile("");
            std::string fname = fvCurve4.getName();
            fname.append(".csv");
            remove(fname.c_str());

        cout << "Passed: Testing Services for connectivity" << endl;                            

       cout <<"Service correctness is tested by underlying utility class"<<endl;
        cout <<"MuscleCurveFunction, and MuscleCurveFunctionFactory"<<endl;

        cout <<"________________________________________________________"<<endl;
        cout <<"          TESTING ForceVelocityCurve              "<<endl;
        cout <<"                    COMPLETED                     "<<endl;
        cout <<"________________________________________________________"<<endl;

}

void testForceVelocityInverseCurve()
{

        //cout <<"**************************************************"<<endl;
        cout <<"________________________________________________________"<<endl;    
        cout <<"1. Testing ForceVelocityInverseCurve"<<endl;       
        cout <<"________________________________________________________"<<endl;

        cout <<"    a. default construction" <<endl;
        ForceVelocityInverseCurve fvCurve1;
        fvCurve1.print("default_ForceVelocityInverseCurve.xml");

        cout <<"    b. serialization & deserialization" <<endl;
        ForceVelocityInverseCurve fvCurve2;
        //change all of the properties to something other than the default
        fvCurve2.setConcentricCurviness(0.5);
        fvCurve2.setConcentricMinSlope(0.05);
        fvCurve2.setEccentricCurviness(0.5);
        fvCurve2.setEccentricMinSlope(0.06);
        fvCurve2.setMaxEccentricVelocityForceMultiplier(2.0);
        fvCurve2.setIsometricMaxSlope(10);


        //These next few lines are just to read the object in, and repopulate
        //fvCurve2 with the properties from the file ... and its a little 
        //awkward to use.

        //cout << "b.*Uncomment, test makeObjectFromFile once in OpenSim"<<endl;        
        
        Object* tmpObj = Object::          
                    makeObjectFromFile("default_ForceVelocityInverseCurve.xml");
        fvCurve2 = *dynamic_cast<ForceVelocityInverseCurve*>(tmpObj);        
        delete tmpObj;
        SimTK_TEST(fvCurve2 == fvCurve1);       
        remove("default_ForceVelocityInverseCurve.xml");
        

        fvCurve2.setConcentricCurviness(0.5);
        fvCurve2.setConcentricMinSlope(0.05);
        fvCurve2.setEccentricCurviness(0.5);
        fvCurve2.setEccentricMinSlope(0.06);
        fvCurve2.setMaxEccentricVelocityForceMultiplier(2.0);
        fvCurve2.setIsometricMaxSlope(10);

        cout <<"    c. assignment operator" <<endl;
        fvCurve2=fvCurve1;
        
        
        SimTK_TEST(fvCurve1==fvCurve2);

        fvCurve2.setConcentricCurviness(0.5);
        fvCurve2.setConcentricMinSlope(0.05);
        fvCurve2.setEccentricCurviness(0.5);
        fvCurve2.setEccentricMinSlope(0.06);
        fvCurve2.setMaxEccentricVelocityForceMultiplier(2.0);
        fvCurve2.setIsometricMaxSlope(10);

        cout <<"    d. copy constructor" <<endl;
        ForceVelocityInverseCurve fvCurve2p5(fvCurve2);
        SimTK_TEST(fvCurve2==fvCurve2p5);

        cout << "*Passed: default construction, limited serialization" << endl;
        cout << "         assignment operator, copy constructor" << endl;

        //cout <<"**************************************************"<<endl;
        cout <<"2. Testing API constructor" << endl;
        ForceVelocityInverseCurve fvCurve3(0.1,5,0.1,1.8,0.1,0.75,"testMuscle");
        double falVal  = fvCurve3.calcValue(1.0);
        double dfalVal = fvCurve3.calcDerivative(1.0,1);
        cout << "Passed: Testing API constructor" << endl;

        //cout <<"**************************************************"<<endl;
        cout <<"3. Testing get/set methods:" << endl;

        fvCurve2.setConcentricCurviness(0.5);
        fvCurve2.setConcentricMinSlope(0.05);
        fvCurve2.setEccentricCurviness(0.5);
        fvCurve2.setEccentricMinSlope(0.06);
        fvCurve2.setMaxEccentricVelocityForceMultiplier(2.0);
        fvCurve2.setIsometricMaxSlope(10);

        SimTK_TEST(fvCurve2.getConcentricCurviness()                    == 0.5);
        SimTK_TEST(fvCurve2.getConcentricMinSlope()                     ==0.05);
        SimTK_TEST(fvCurve2.getEccentricCurviness()                     == 0.5);
        SimTK_TEST(fvCurve2.getEccentricMinSlope()                      ==0.06);
        SimTK_TEST(fvCurve2.getMaxEccentricVelocityForceMultiplier()    == 2.0);
        SimTK_TEST(fvCurve2.getIsometricMaxSlope()                      ==  10);

        cout << "Passed: Testing get/set methods" << endl;

//====================================================================
        double p1 = 0.1;
        double p2 = 5;
        double p3 = 0.1;
        double p4 = 1.8;
        double p5 = 0.1;
        double p6 = 0.75;

        printf("4. Testing default property values:"
               " \n\t%f,\n\t%f,\n\t%f,\n\t%f,\n\t%f,\n\t%f\n"
            ,p1,p2,p3,p4,p5,p6);      
        ForceVelocityInverseCurve fvCurve4;
        fvCurve4.setName("fvInvCurve");

        SimTK_TEST(fvCurve4.getConcentricMinSlope()     == p1);
        SimTK_TEST(fvCurve4.getIsometricMaxSlope()      == p2);
        SimTK_TEST(fvCurve4.getEccentricMinSlope()      == p3);
        SimTK_TEST(fvCurve4.getMaxEccentricVelocityForceMultiplier() == p4);
        SimTK_TEST(fvCurve4.getConcentricCurviness()    == p5);
        SimTK_TEST(fvCurve4.getEccentricCurviness()     == p6);
        cout << "Passed" << endl;
//====================================================================
        cout <<"5. Testing serivices for connectivity:" << endl;     


        cout <<"    a. calcValue" << endl;
            double tol = sqrt(SimTK::Eps);
            double value = fvCurve4.calcValue(1.0);
            SimTK_TEST_EQ_TOL(value, 0, tol);
        cout <<"    b. calcDerivative" << endl;
            double dvalue= fvCurve4.calcDerivative(1.0,1);
            SimTK_TEST_EQ_TOL(dvalue, 1.0/5.0, tol);
            dvalue= fvCurve4.calcDerivative(1.0,2);
            SimTK_TEST_EQ_TOL(dvalue, 0, tol);

        cout <<"    c. getCurveDomain" << endl;
            SimTK::Vec2 tmp = fvCurve4.getCurveDomain();
            SimTK_TEST(tmp(0) == 0 &&
                       tmp(1) == p4);

        cout <<"    d. printMuscleCurveToCSVFile" << endl;
            fvCurve4.setConcentricCurviness(0.5);
            fvCurve4.setEccentricCurviness(1.0);
            fvCurve4.printMuscleCurveToCSVFile("");
            std::string fname = fvCurve4.getName();
            fname.append(".csv");
            remove(fname.c_str());

        cout << "Passed: Testing Services for connectivity" << endl;                            

        //cout <<"**************************************************"<<endl;
       cout <<"Service correctness is tested by underlying utility class"<<endl;
       cout <<"MuscleCurveFunction, and MuscleCurveFunctionFactory"<<endl;
        //cout <<"**************************************************"<<endl;

        //cout <<"**************************************************"<<endl;
        cout <<"________________________________________________________"<<endl;
        cout <<"          TESTING ForceVelocityInverseCurve             "<<endl;
        cout <<"                    COMPLETED                     "<<endl;
        cout <<"________________________________________________________"<<endl;
        //cout <<"**************************************************"<<endl;

}


void testTendonForceLengthCurve()
{

        cout <<"________________________________________________________"<<endl;    
        cout <<"1. Testing TendonForceLengthCurve"<<endl;       
        cout <<"________________________________________________________"<<endl;

        cout <<"    a. default construction" <<endl;
        TendonForceLengthCurve fseCurve1;
        fseCurve1.print("default_TendonForceLengthCurve.xml");

        cout <<"    b. serialization & deserialization" <<endl;
        TendonForceLengthCurve fseCurve2;
        //change all of the properties to something other than the default
        fseCurve2.setStrainAtOneNormForce(0.10);
        fseCurve2.setStiffnessAtOneNormForce(50.0);
        fseCurve2.setCurviness(0.8);


        //These next few lines are just to read the object in, and repopulate
        //fvCurve2 with the properties from the file ... and its a little 
        //awkward to use.

        //cout << "b.*Uncomment, test makeObjectFromFile once in OpenSim"<<endl;        
        
        Object* tmpObj = Object::          
                    makeObjectFromFile("default_TendonForceLengthCurve.xml");
        fseCurve2 = *dynamic_cast<TendonForceLengthCurve*>(tmpObj);        
        delete tmpObj;
        SimTK_TEST(fseCurve2 == fseCurve1);       
        remove("default_TendonForceLengthCurve.xml");
        

        fseCurve2.setStrainAtOneNormForce(0.10);
        fseCurve2.setStiffnessAtOneNormForce(50.0);
        fseCurve2.setCurviness(0.8);

        cout <<"    c. assignment operator" <<endl;
        fseCurve2=fseCurve1;
                
        SimTK_TEST(fseCurve1==fseCurve2);

        fseCurve2.setStrainAtOneNormForce(0.10);
        fseCurve2.setStiffnessAtOneNormForce(50.0);
        fseCurve2.setCurviness(0.8);

        cout <<"    d. copy constructor" <<endl;
        TendonForceLengthCurve fseCurve2p5(fseCurve2);
        SimTK_TEST(fseCurve2==fseCurve2p5);

        cout << "*Passed: default construction, limited serialization" << endl;
        cout << "         assignment operator, copy constructor" << endl;

        //====================================================================
        cout <<"2. Testing API constructor" << endl;
        TendonForceLengthCurve fseCurve3(0.10,50,0.75,"testMuscle");
        double falVal  = fseCurve3.calcValue(0.02);
        double dfalVal = fseCurve3.calcDerivative(0.02,1);
        cout << "Passed: Testing API constructor" << endl;

        //====================================================================
        cout <<"3. Testing get/set methods:" << endl;

        fseCurve2.setStrainAtOneNormForce(0.10);
        fseCurve2.setStiffnessAtOneNormForce(50.0);
        fseCurve2.setCurviness(0.8);

        SimTK_TEST(fseCurve2.getStrainAtOneNormForce()      == 0.10);
        SimTK_TEST(fseCurve2.getStiffnessAtOneNormForce()   == 50.0);
        SimTK_TEST(fseCurve2.getCurviness()                 == 0.80);
        
        cout << "Passed: Testing get/set methods" << endl;

        //====================================================================
        double p1 = 0.04;
        double p2 = 42;
        double p3 = 0.75;

        printf("4. Testing default property values: \n\t%f,\n\t%f,\n\t%f\n"
            ,p1,p2,p3);
        TendonForceLengthCurve fseCurve4;
            SimTK_TEST(fseCurve4.getStrainAtOneNormForce()      == p1);
            SimTK_TEST(fseCurve4.getStiffnessAtOneNormForce()   == p2);
            SimTK_TEST(fseCurve4.getCurviness()                 == p3);
        cout << "Passed" << endl;
        //====================================================================
        cout <<"5. Testing default curve values at end points,"
               " and services" << endl;
        fseCurve4.setName("fseCurve");

        cout <<"    a. calcValue" << endl;
            double l0 = 1;
            double l1 = l0 + p1;
            double dydx = p2;

            double tol = sqrt(SimTK::Eps);

            double value = fseCurve4.calcValue(l0);
                SimTK_TEST_EQ_TOL(value, 0, tol);
            value = fseCurve4.calcValue(l1);
                SimTK_TEST_EQ_TOL(value, 1, tol);
        cout <<"    b. calcDerivative" << endl;
            double dvalue= fseCurve4.calcDerivative(l0,1);
                SimTK_TEST_EQ_TOL(dvalue, 0, tol);
            dvalue= fseCurve4.calcDerivative(l1,1);
                SimTK_TEST_EQ_TOL(dvalue, dydx, tol);
            //look at the second derivative
            dvalue= fseCurve4.calcDerivative(l0,2);
                SimTK_TEST_EQ_TOL(dvalue, 0, tol);
            dvalue= fseCurve4.calcDerivative(l1,2);
                SimTK_TEST_EQ_TOL(dvalue, 0, tol);

        cout <<"    c. getCurveDomain" << endl;
            SimTK::Vec2 tmp = fseCurve4.getCurveDomain();
            SimTK_TEST(tmp(0) == l0 &&
                       tmp(1) == l1);

        cout <<"    d. printMuscleCurveToCSVFile" << endl;            

            fseCurve4.printMuscleCurveToCSVFile("");
            std::string fname = fseCurve4.getName();
            fname.append(".csv");
            remove(fname.c_str());

       cout << "Passed: Testing Services for connectivity" << endl;                            

       cout <<"Service correctness is tested by underlying utility class"<<endl;
       cout <<"MuscleCurveFunction, and MuscleCurveFunctionFactory"<<endl;

       cout <<"________________________________________________________"<<endl;
       cout <<"          TESTING TendonForceLengthCurve             "<<endl;
       cout <<"                    COMPLETED                     "<<endl;
       cout <<"________________________________________________________"<<endl;

}

void testFiberForceLengthCurve()
{

        cout <<"________________________________________________________"<<endl;    
        cout <<"1. Testing FiberForceLengthCurve"<<endl;       
        cout <<"________________________________________________________"<<endl;

        cout <<"    a. default construction" <<endl;
        FiberForceLengthCurve fpeCurve1;
        fpeCurve1.print("default_FiberForceLengthCurve.xml");

        cout <<"    b. serialization & deserialization" <<endl;
        FiberForceLengthCurve fpeCurve2;
        //change all of the properties to something other than the default
        fpeCurve2.setStrainAtOneNormForce(0.80);
        fpeCurve2.setStiffnessAtOneNormForce(10.0);
        fpeCurve2.setCurviness(0.8);


        //These next few lines are just to read the object in, and repopulate
        //fvCurve2 with the properties from the file ... and its a little 
        //awkward to use.

        //cout << "b.*Uncomment, test makeObjectFromFile once in OpenSim"<<endl;        
        
        Object* tmpObj = Object::          
                    makeObjectFromFile("default_FiberForceLengthCurve.xml");
        fpeCurve2 = *dynamic_cast<FiberForceLengthCurve*>(tmpObj);        
        delete tmpObj;
        SimTK_TEST(fpeCurve2 == fpeCurve1);       
        remove("default_FiberForceLengthCurve.xml");
        

        fpeCurve2.setStrainAtOneNormForce(0.80);
        fpeCurve2.setStiffnessAtOneNormForce(10.0);
        fpeCurve2.setCurviness(0.8);

        cout <<"    c. assignment operator" <<endl;
        fpeCurve2=fpeCurve1;
                
        SimTK_TEST(fpeCurve1==fpeCurve2);

        fpeCurve2.setStrainAtOneNormForce(0.80);
        fpeCurve2.setStiffnessAtOneNormForce(10.0);
        fpeCurve2.setCurviness(0.8);

        cout <<"    d. copy constructor" <<endl;
        FiberForceLengthCurve fpeCurve2p5(fpeCurve2);
        SimTK_TEST(fpeCurve2==fpeCurve2p5);

        cout << "*Passed: default construction, limited serialization" << endl;
        cout << "         assignment operator, copy constructor" << endl;

        //====================================================================
        cout <<"2. Testing API constructor" << endl;
        FiberForceLengthCurve fpeCurve3(0.10,50,0.75,"testMuscle");
        double falVal  = fpeCurve3.calcValue(0.02);
        double dfalVal = fpeCurve3.calcDerivative(0.02,1);
        cout << "Passed: Testing API constructor" << endl;

        //====================================================================
        cout <<"3. Testing get/set methods:" << endl;

        fpeCurve2.setStrainAtOneNormForce(0.80);
        fpeCurve2.setStiffnessAtOneNormForce(10.0);
        fpeCurve2.setCurviness(0.8);

        SimTK_TEST(fpeCurve2.getStrainAtOneNormForce()      == 0.80);
        SimTK_TEST(fpeCurve2.getStiffnessAtOneNormForce()   == 10.0);
        SimTK_TEST(fpeCurve2.getCurviness()                 == 0.80);
        
        cout << "Passed: Testing get/set methods" << endl;

        //====================================================================
        double p1 = 0.6;
        double p2 = 8.4;
        double p3 = 0.65;

        printf("4. Testing default property values: \n\t%f,\n\t%f,\n\t%f\n",
            p1,p2,p3);      
        FiberForceLengthCurve fpeCurve4;
        fpeCurve4.setName("fpeCurve");

        SimTK_TEST(fpeCurve4.getStrainAtOneNormForce()      == p1);
        SimTK_TEST(fpeCurve4.getStiffnessAtOneNormForce()   == p2);
        SimTK_TEST(fpeCurve4.getCurviness()                 == p3);

        fpeCurve4.setName("fpeCurve");
        cout << "Passed" << endl;
        //====================================================================
        cout <<"5. Testing default curve values at end points,"
               " and services" << endl;

        //cout <<"    a. calcValue" << endl;
            double l0 = 1.0;
            double l1 = l0+p1;
            double dydx = p2;

        cout <<"    a. calcValue" << endl;
            double tol = sqrt(SimTK::Eps);
            double value = fpeCurve4.calcValue(l0);
            SimTK_TEST_EQ_TOL(value, 0, tol);
            value = fpeCurve4.calcValue(l1);
            SimTK_TEST_EQ_TOL(value, 1, tol);
        cout <<"    b. calcDerivative" << endl;
            double dvalue= fpeCurve4.calcDerivative(l0,1);
            SimTK_TEST_EQ_TOL(dvalue, 0, tol);
            dvalue= fpeCurve4.calcDerivative(l1,1);
            SimTK_TEST_EQ_TOL(dvalue, dydx, tol);
            
            //look at the second derivative
            dvalue= fpeCurve4.calcDerivative(l0,2);
            SimTK_TEST_EQ_TOL(dvalue, 0, tol);
            dvalue= fpeCurve4.calcDerivative(l1,2);
            SimTK_TEST_EQ_TOL(dvalue, 0, tol);

        cout <<"    c. getCurveDomain" << endl;
            SimTK::Vec2 tmp = fpeCurve4.getCurveDomain();
            SimTK_TEST(tmp(0) == l0 &&
                       tmp(1) == l1);

        cout <<"    d. printMuscleCurveToCSVFile" << endl;            

            fpeCurve4.printMuscleCurveToCSVFile("");
            std::string fname = fpeCurve4.getName();
            fname.append(".csv");
            remove(fname.c_str());

        cout << "Passed: Testing Services for connectivity" << endl;                            

       cout <<"Service correctness is tested by underlying utility class"<<endl;
       cout <<"MuscleCurveFunction, and MuscleCurveFunctionFactory"<<endl;

        cout <<"________________________________________________________"<<endl;
        cout <<"          TESTING FiberForceLengthCurve             "<<endl;
        cout <<"                    COMPLETED                     "<<endl;
        cout <<"________________________________________________________"<<endl;

}

void testFiberCompressiveForceLengthCurve()
{

        cout <<"________________________________________________________"<<endl;    
        cout <<"1. Testing FiberCompressiveForceLengthCurve"<<endl;       
        cout <<"________________________________________________________"<<endl;

        cout <<"    a. default construction" <<endl;
        FiberCompressiveForceLengthCurve fceCurve1;
        fceCurve1.print("default_FiberCompressiveForceLengthCurve.xml");

        cout <<"    b. serialization & deserialization" <<endl;
        FiberCompressiveForceLengthCurve fceCurve2;
        //change all of the properties to something other than the default
        fceCurve2.setNormLengthAtZeroForce(0.7);
        fceCurve2.setStiffnessAtZeroLength(-10.0);
        fceCurve2.setCurviness(0.1);


        //These next few lines are just to read the object in, and repopulate
        //fvCurve2 with the properties from the file ... and its a little 
        //awkward to use.

        //cout << "b.*Uncomment, test makeObjectFromFile once in OpenSim"<<endl;        
        
        Object* tmpObj = Object::          
             makeObjectFromFile("default_FiberCompressiveForceLengthCurve.xml");
        fceCurve2 = *dynamic_cast<FiberCompressiveForceLengthCurve*>(tmpObj);        
        delete tmpObj;
        SimTK_TEST(fceCurve2 == fceCurve1);       
        remove("default_FiberCompressiveForceLengthCurve.xml");
        

        //change all of the properties to something other than the default
        fceCurve2.setNormLengthAtZeroForce(0.7);
        fceCurve2.setStiffnessAtZeroLength(-10.0);
        fceCurve2.setCurviness(0.1);

        cout <<"    c. assignment operator" <<endl;
        fceCurve2=fceCurve1;
                
        SimTK_TEST(fceCurve1==fceCurve2);

        //change all of the properties to something other than the default
        fceCurve2.setNormLengthAtZeroForce(0.7);
        fceCurve2.setStiffnessAtZeroLength(-10.0);
        fceCurve2.setCurviness(0.1);

        cout <<"    d. copy constructor" <<endl;
        FiberCompressiveForceLengthCurve fceCurve2p5(fceCurve2);
        SimTK_TEST(fceCurve2==fceCurve2p5);

        cout << "*Passed: default construction, limited serialization" << endl;
        cout << "         assignment operator, copy constructor" << endl;

        //====================================================================
        cout <<"2. Testing API constructor" << endl;
        FiberCompressiveForceLengthCurve fceCurve3(0.80,-10,0.1,"testMuscle");
        double fceVal  = fceCurve3.calcValue(0.02);
        double dfceVal = fceCurve3.calcDerivative(0.02,1);
        cout << "Passed: Testing API constructor" << endl;

        //====================================================================
        cout <<"3. Testing get/set methods:" << endl;

        //change all of the properties to something other than the default
        fceCurve2.setNormLengthAtZeroForce(0.7);
        fceCurve2.setStiffnessAtZeroLength(-10.0);
        fceCurve2.setCurviness(0.1);

        SimTK_TEST(fceCurve2.getNormLengthAtZeroForce()     == 0.70);
        SimTK_TEST(fceCurve2.getStiffnessAtZeroLength()     == -10.0);
        SimTK_TEST(fceCurve2.getCurviness()                 == 0.10);
        
        cout << "Passed: Testing get/set methods" << endl;

        //====================================================================
        double p1 = 0.58564173314080115;
        double p2 = -8.0;
        double p3 = 0.1;

        printf("4a. Testing default property values: \n\t%f,\n\t%f,\n\t%f\n",
            p1,p2,p3);      
        FiberCompressiveForceLengthCurve fceCurve4;
        fceCurve4.setName("fceCurve");

        SimTK_TEST(fceCurve4.getNormLengthAtZeroForce()     == p1);
        SimTK_TEST(fceCurve4.getStiffnessAtZeroLength()     == p2);
        SimTK_TEST(fceCurve4.getCurviness()                 == p3);
        cout << "Passed" << endl;

        //solve for the paramters that will yield a value of 0.05 at a 
        //fiber length of  0.5

        double len0 = fceCurve4.getNormLengthAtZeroForce();
        double desLen = 0.5;
        double desVal = 0.05;
        double tol = 1e-12;
        /*
        double step = 0.1;
        int maxIter = 100;
        int iter = 0;
       
        double err      = abs(desVal - fceCurve4.calcValue(desLen));
        double errLeft  = 1;
        double errRight = 1;

        while(abs(err) > tol && iter < maxIter){
            fceCurve4.setNormLengthAtZeroForce(len0 - step);
            errLeft = desVal - fceCurve4.calcValue(desLen);
            
            fceCurve4.setNormLengthAtZeroForce(len0 + step);
            errRight = desVal - fceCurve4.calcValue(desLen);

            if(abs(errLeft) < abs(err) || abs(errRight) < abs(err)){
                if(abs(errLeft) < abs(errRight)){
                    err = abs(errLeft);
                    len0 = len0-step;
                }else{
                    err = abs(errRight);
                    len0 = len0+step;
                }
            }
            step = step/2;
            iter = iter+1;
        }
        */
        //printf("NormLengthAtZeroForce %f will yield %f at %f\n",
        //        len0,desVal,desLen);
        
        tol = sqrt(SimTK::Eps);        
        printf("4b. Checking that at %f the curve value is %f\n",desLen,desVal);
        SimTK_TEST_EQ_TOL(fceCurve4.calcValue(desLen), desVal, tol);
        cout << "Passed" << endl;

        fceCurve4.setName("fceCurve");

        //====================================================================
        cout <<"5. Testing default curve values at end points,"
               " and services" << endl;

        //cout <<"    a. calcValue" << endl;
            double l0 = p1;
            double l1 = 0;
            double dydx = p2;

        cout <<"    a. calcValue" << endl;
            tol = sqrt(SimTK::Eps);
            double value = fceCurve4.calcValue(l0);
            SimTK_TEST_EQ_TOL(value, 0, tol);
            value = fceCurve4.calcValue(l1);
            SimTK_TEST_EQ_TOL(value, 1, tol);
        cout <<"    b. calcDerivative" << endl;
            double dvalue= fceCurve4.calcDerivative(l0,1);
            SimTK_TEST_EQ_TOL(dvalue, 0, tol);
            dvalue= fceCurve4.calcDerivative(l1,1);
            SimTK_TEST_EQ_TOL(dvalue, dydx, tol);
            
            //look at the second derivative
            dvalue= fceCurve4.calcDerivative(l0,2);
            SimTK_TEST_EQ_TOL(dvalue, 0, tol);
            dvalue= fceCurve4.calcDerivative(l1,2);
            SimTK_TEST_EQ_TOL(dvalue, 0, tol);

        cout <<"    c. getCurveDomain" << endl;
            SimTK::Vec2 tmp = fceCurve4.getCurveDomain();
            SimTK_TEST(tmp(0) == l1 &&
                       tmp(1) == l0);




        cout <<"    d. printMuscleCurveToCSVFile" << endl;            

            fceCurve4.printMuscleCurveToCSVFile("");
            std::string fname = fceCurve4.getName();
            fname.append(".csv");
            remove(fname.c_str());

        cout << "Passed: Testing Services for connectivity" << endl;                            

       cout <<"Service correctness is tested by underlying utility class"<<endl;
       cout <<"MuscleCurveFunction, and MuscleCurveFunctionFactory"<<endl;

        cout <<"________________________________________________________"<<endl;
        cout <<"          TESTING FiberCompressiveForceLengthCurve      "<<endl;
        cout <<"                    COMPLETED                     "<<endl;
        cout <<"________________________________________________________"<<endl;

}

void testFiberCompressiveForceCosPennationCurve()
{

        cout <<"________________________________________________________"<<endl;    
        cout <<"1. Testing FiberCompressiveForceCosPennationCurve"<<endl;       
        cout <<"________________________________________________________"<<endl;

        cout <<"    a. default construction" <<endl;
        FiberCompressiveForceCosPennationCurve fcpCurve1;
        fcpCurve1.print("default_FiberCompressiveForceCosPennationCurve.xml");

        cout <<"    b. serialization & deserialization" <<endl;
        FiberCompressiveForceCosPennationCurve fcpCurve2;
        //change all of the properties to something other than the default
        fcpCurve2.setEngagementAngleInDegrees(85.0);
        fcpCurve2.setStiffnessAtPerpendicular(-20.0);
        fcpCurve2.setCurviness(0.5);


        //These next few lines are just to read the object in, and repopulate
        //fvCurve2 with the properties from the file ... and its a little 
        //awkward to use.

        cout << "b.*Uncomment, test makeObjectFromFile once in OpenSim"<<endl;        
        
        Object* tmpObj = Object::          
                    makeObjectFromFile(
                    "default_FiberCompressiveForceCosPennationCurve.xml");
        fcpCurve2=*dynamic_cast<FiberCompressiveForceCosPennationCurve*>(tmpObj);        
        delete tmpObj;
        SimTK_TEST(fcpCurve2 == fcpCurve1);       
        remove("default_FiberCompressiveForceCosPennationCurve.xml");
        

        //change all of the properties to something other than the default
        fcpCurve2.setEngagementAngleInDegrees(85.0);
        fcpCurve2.setStiffnessAtPerpendicular(-20.0);
        fcpCurve2.setCurviness(0.5);

        cout <<"    c. assignment operator" <<endl;
        fcpCurve2=fcpCurve1;
                
        SimTK_TEST(fcpCurve1==fcpCurve2);

        //change all of the properties to something other than the default
        fcpCurve2.setEngagementAngleInDegrees(85.0);
        fcpCurve2.setStiffnessAtPerpendicular(-20.0);
        fcpCurve2.setCurviness(0.5);

        cout <<"    d. copy constructor" <<endl;
        FiberCompressiveForceCosPennationCurve fcpCurve2p5(fcpCurve2);
        SimTK_TEST(fcpCurve2==fcpCurve2p5);

        cout << "*Passed: default construction, limited serialization" << endl;
        cout << "         assignment operator, copy constructor" << endl;

        //====================================================================
        cout <<"2. Testing API constructor(s)" << endl;
        FiberCompressiveForceCosPennationCurve 
            fcpCurve3(78,-10,0.1,"testMuscle");
        
        double cosAngle = cos(85*SimTK::Pi/180);

        double fcpVal  = fcpCurve3.calcValue(cosAngle);
        double dfcpVal = fcpCurve3.calcDerivative(cosAngle,1);
        
        FiberCompressiveForceCosPennationCurve fcpCurve3a(78,"testMuscle");
        fcpVal  = fcpCurve3a.calcValue(cosAngle);
        dfcpVal = fcpCurve3a.calcDerivative(cosAngle,1);

        FiberCompressiveForceCosPennationCurve fcpCurve3b;
        fcpVal  = fcpCurve3b.calcValue(cosAngle);
        dfcpVal = fcpCurve3b.calcDerivative(cosAngle,1);

        cout << "Passed: Testing API constructor" << endl;

        //====================================================================
        cout <<"3. Testing get/set methods:" << endl;

        //change all of the properties to something other than the default
        fcpCurve2.setEngagementAngleInDegrees(85.0);
        fcpCurve2.setStiffnessAtPerpendicular(-20.0);
        fcpCurve2.setCurviness(0.5);

        SimTK_TEST(fcpCurve2.getEngagementAngleInDegrees()  ==  85.0);
        SimTK_TEST(fcpCurve2.getStiffnessAtPerpendicular()  == -20.0);
        SimTK_TEST(fcpCurve2.getCurviness()                 == 0.50);
       
        cout << "Passed: Testing get/set methods" << endl;

        //====================================================================
        double p1 = 80;
        double p2 = -2.0/cos(p1*SimTK::Pi/180.0);
        double p3 = 0.1;

        printf("4a. Testing default property values: \n\t%f,\n\t%f (computed)"
            ",\n\t%f (computed)\n",
            p1,p2,p3);      
        FiberCompressiveForceCosPennationCurve fcpCurve4;
        fcpCurve4.setName("fceCurve");

        SimTK_TEST(fcpCurve4.getEngagementAngleInDegrees()      == p1);
        SimTK_TEST(fcpCurve4.getStiffnessAtPerpendicularInUse() == p2);
        SimTK_TEST(fcpCurve4.getCurvinessInUse()                == p3);
        cout << "Passed" << endl;

        fcpCurve4.setName("fcpCurve");

        //====================================================================
        cout <<"5. Testing default curve values at end points,"
               " and services" << endl;

        //cout <<"    a. calcValue" << endl;
            double l0 = cos(p1*SimTK::Pi/180);
            double l1 = cos(90*SimTK::Pi/180);
            double dydx = p2;

        cout <<"    a. calcValue" << endl;
            double tol = sqrt(SimTK::Eps);
            double value = fcpCurve4.calcValue(l0);
            SimTK_TEST_EQ_TOL(value, 0, tol);
            value = fcpCurve4.calcValue(l1);
            SimTK_TEST_EQ_TOL(value, 1, tol);

        cout <<"    b. calcDerivative" << endl;
            double dvalue= fcpCurve4.calcDerivative(l0,1);
            SimTK_TEST_EQ_TOL(dvalue, 0, tol);
            dvalue= fcpCurve4.calcDerivative(l1,1);
            SimTK_TEST_EQ_TOL(dvalue, dydx, tol);
            
            //look at the second derivative
            dvalue= fcpCurve4.calcDerivative(l0,2);
            SimTK_TEST_EQ_TOL(dvalue, 0, tol);
            dvalue= fcpCurve4.calcDerivative(l1,2);
            SimTK_TEST_EQ_TOL(dvalue, 0, tol);

        cout <<"    c. getCurveDomain" << endl;
            SimTK::Vec2 tmp = fcpCurve4.getCurveDomain();
            SimTK_TEST_EQ_TOL(tmp(0), l1, tol); 
            SimTK_TEST_EQ_TOL(tmp(1), l0, tol); 




        cout <<"    d. printMuscleCurveToCSVFile" << endl;            

            fcpCurve4.printMuscleCurveToCSVFile("");
            std::string fname = fcpCurve4.getName();
            fname.append(".csv");
            remove(fname.c_str());

        cout << "Passed: Testing Services for connectivity" << endl;                            

       cout <<"Service correctness is tested by underlying utility class"<<endl;
       cout <<"MuscleCurveFunction, and MuscleCurveFunctionFactory"<<endl;

        cout <<"________________________________________________________"<<endl;
        cout <<"          TESTING FiberCompressiveForceCosPennationCurve"<<endl;
        cout <<"                    COMPLETED                     "<<endl;
        cout <<"________________________________________________________"<<endl;

}
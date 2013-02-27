#include <iostream>
#include "Vector.h"
#include "math.h"
#include "FormFactorData.h"
//#include "LatticePlane.h"
#include "XRay.h"
#include "AbsorbCoeffData.h"
#include "PowderDiffraction.h"
//#include <random> need this on windows?
#include <float.h>
#include <time.h>
//#include <boost/random/linear_congruential.hpp>
//#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include "CCD.h"
#include "FileReading.h"

using namespace std;

#define FLU_ENABLE

#ifndef base_generator_type
typedef boost::mt11213b base_generator_type;
#endif


//http://www.boost.org/doc/libs/1_46_1/libs/random/example/random_demo.cpp
//http://www.boost.org/doc/libs/1_52_0/doc/html/boost_random/reference.html

int main()
{
    ifstream datafile("InputScript.txt");
    std::map<std::string,std::string> InputData;
    AddToMapFromFile(datafile, InputData);
    datafile.close();
    //exit(0);

    Vector InputCCDOrigin(110,-1,50); //origin of CCD
    Vector InputCCDNormal(0,0,1); //direction that CCD points in.
    double InputCCDAngle = 0;

    double InputCCDXMin = 0.0;
    double InputCCDXMax = 5.0;

    double InputCCDYMin = 0.0;
    double InputCCDYMax = 2.0;

    VectorFromMap("CCDOrigin",InputData,InputCCDOrigin);
    VectorFromMap("CCDNormal",InputData,InputCCDNormal);
    DoubleFromMap("CCDAngle", InputData,InputCCDAngle);
    DoubleFromMap("CCDXMin", InputData, InputCCDXMin);
    DoubleFromMap("CCDYMin", InputData, InputCCDYMin);
    DoubleFromMap("CCDXMax", InputData, InputCCDXMax);
    DoubleFromMap("CCDYMax", InputData, InputCCDYMax);

    InputCCDOrigin.Print();
    InputCCDNormal.Print();


    CCD CCDCamera(InputCCDOrigin, InputCCDNormal, InputCCDAngle,
                  0.05, 0.05,
                  InputCCDXMin, InputCCDXMax,
                  InputCCDYMin, InputCCDYMax);

    //Vector CCDCorner1, CCDCorner2, CCDCorner3, CCDCorner4;

    Vector CCDCorners[4];

    CCDCamera.GenerateCCDCorners( CCDCorners[0], CCDCorners[1], CCDCorners[2], CCDCorners[3]);

    cout << "CCDBounds:" << endl;
    CCDCorners[0].Print();
    CCDCorners[1].Print();
    CCDCorners[2].Print();
    CCDCorners[3].Print();

    Vector CrystalCorner1(0,-0.1,0), CrystalCorner2(0,0.1,0), CrystalCorner3(0.1,-0.1,0), CrystalCorner4(0.1,0.1,0);

    Vector Corners[4];
    Corners[0] = CrystalCorner1;
    Corners[1] = CrystalCorner2;
    Corners[2] = CrystalCorner3;
    Corners[3] = CrystalCorner4;

    double DirectionMinX = 0.0, DirectionMinY = 0.0, DirectionMaxX = 0.0, DirectionMaxY = 0.0;
    double DirectionMinCosTheta = 0.0, DirectionMaxCosTheta = 0.0, DirectionMinPhi = 0.0, DirectionMaxPhi = 0.0;


    bool bFirstIteration = true;

    //loop through each corner of the crystal, then loop through each of the CCD corners to find the min/max for the solid angle rejection parameters.
    for(int iCrystalCorner = 0; iCrystalCorner<4; iCrystalCorner++)
    {
        for( int iCCDCorner = 0; iCCDCorner<4; iCCDCorner++)
        {
            Vector Direction = CCDCorners[iCCDCorner] - Corners[iCrystalCorner];
            Direction = Direction.Normalized();

            double CosTheta = cos(Direction.GetTheta());
            double Phi = Direction.GetPhi();

            if(bFirstIteration)
            {
                DirectionMinX = Direction.x;
                DirectionMaxX = Direction.x;
                DirectionMinY = Direction.y;
                DirectionMaxY = Direction.y;

                DirectionMinCosTheta = CosTheta;
                DirectionMaxCosTheta = CosTheta;
                DirectionMinPhi = Phi;
                DirectionMaxPhi = Phi;

                bFirstIteration = false;
                continue;
            }

            if(Direction.x < DirectionMinX)
            {
                DirectionMinX = Direction.x;
            }
            if(Direction.x > DirectionMaxX)
            {
                DirectionMaxX = Direction.x;
            }
            if(Direction.y < DirectionMinY)
            {
                DirectionMinY = Direction.y;
            }
            if(Direction.y > DirectionMaxY)
            {
                DirectionMaxY = Direction.y;
            }

            if(Phi < DirectionMinPhi)
            {
                DirectionMinPhi = Phi;
            }
            if(Phi > DirectionMaxPhi)
            {
                DirectionMaxPhi = Phi;
            }
            if(CosTheta < DirectionMinCosTheta)
            {
                DirectionMinCosTheta = CosTheta;
            }
            if(CosTheta > DirectionMaxCosTheta)
            {
                DirectionMaxCosTheta = CosTheta;
            }
        }
    }

    /*cout << "X:\t" << DirectionMinX << "\t" << DirectionMaxX << endl;
    cout << "Y:\t" << DirectionMinY << "\t" << DirectionMaxY << endl;*/

    double SolidAngleDelta = 0.005; //% tolerance for solid angle rejection. (how much we add/take away to the x&y limits)

    DoubleFromMap("SolidAngleDelta", InputData, SolidAngleDelta);

    DirectionMinX -= DirectionMinX*SolidAngleDelta;
    DirectionMinY -= DirectionMinY*SolidAngleDelta;
    DirectionMaxX += DirectionMaxX*SolidAngleDelta;
    DirectionMaxY += DirectionMaxY*SolidAngleDelta;

    DirectionMinCosTheta-= DirectionMinCosTheta*SolidAngleDelta;
    DirectionMinPhi -= DirectionMinPhi*SolidAngleDelta;
    DirectionMaxCosTheta += DirectionMaxCosTheta*SolidAngleDelta;
    DirectionMaxPhi += DirectionMaxPhi*SolidAngleDelta;


    cout << "X:\t" << DirectionMinX << "\t" << DirectionMaxX << endl;
    cout << "Y:\t" << DirectionMinY << "\t" << DirectionMaxY << endl;

    cout << "CosTheta:\t" << DirectionMinCosTheta << "\t" << DirectionMaxCosTheta << endl;
    cout << "Phi:\t" << DirectionMinPhi << "\t" << DirectionMaxPhi << endl;


    base_generator_type generator(42u);
    boost::uniform_real<> uni_dist(0,1);
    boost::normal_distribution<> normal_dist(0,1);
    boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
    boost::variate_generator<base_generator_type&, boost::normal_distribution<> > normal(generator, normal_dist);

#ifndef FLU_ENABLE
    cout << "Fluorescence Disabled!" << endl;
#endif

    double LatticeConst;

    DoubleFromMap("LatticeConstant", InputData, LatticeConst);

    PowderDiffraction PD( LatticeConst, &uni ); //Tantalum a0 = 3.31A

    std::string DiffractionData;
    StringFromMap("DiffractionData", InputData, DiffractionData);
    PD.LoadData(DiffractionData.c_str());

    AbsorbCoeffData TaMuData( 1.0f, 15.0f, 2000);

    std::string AbsorptionData;
    StringFromMap("AbsorptionData", InputData, AbsorptionData);
    TaMuData.LoadData(AbsorptionData.c_str());

    double Thickness = 500000;//try 50 microns //500000.0f; //50 microns in A

    DoubleFromMap("CrystalThickness", InputData, Thickness);

    //simple detector at X/Y plane
    ofstream DiffractResults( "DiffractResults.txt" );
    ofstream FluoResults( "FluoResults.txt" );

    //full ray information
    ofstream AdvDiffractResults( "AdvDiffractResults.txt");
    ofstream AdvFluoResults( "AdvFluoResults.txt" );


    double SourceZ = 5.0;
    double CCDZ = 50.0;

    DoubleFromMap("SourceZ", InputData, SourceZ);
    DoubleFromMap("CCDZ", InputData, CCDZ);

    Vector Source( -3.0f, 0.0f, SourceZ);

    double MinE = 3.0;//4.23;
    double MaxE = 9.0;//4.28;
    double DeltaE = 0.0005;//0.0005f;

    DoubleFromMap("MinEnergy", InputData, MinE);
    DoubleFromMap("MaxEnergy", InputData, MaxE);
    DoubleFromMap("DeltaE", InputData, DeltaE);

    int nEPoints = (MaxE - MinE)/DeltaE;
    //float Originalx = Source.x;

    long unsigned int StartTime = time(NULL);

    int nDiffracted = 0, nFluoresced = 0;

    for(float Sx = 0.0f; Sx <= 0.0f; Sx += 0.5f)
    {
        //Source.x = Sx+Originalx;
        float ProgressCounter = 0.0f;

        for( int EnergyTick = 0; EnergyTick < nEPoints; EnergyTick++)
        {
            double Energy = MinE + DeltaE*EnergyTick;
            float ProbScatter = PD.GetModifiedScatterProb(Energy);
            float AbsorbCoeff = TaMuData.GetAbsorbCoeffDataPoint( EnergyToWavelength(Energy) );
            for( double x = 3; x <= 3.1; x += 0.01)
            {
                for( double y = -0.1; y <= 0.1; y += 0.01)
                {
                    //construct a vector that goes down SourceZ units and across a variable amount
                    Vector Direction( x, y, -1.0*SourceZ);
                    Direction = Direction.Normalized();

                    float PathLength = fabs(Thickness/Direction.z); //length xray takes through crystal
                    float ProbAbsorb = 1.0f - exp( -1.0f * AbsorbCoeff * PathLength);

                    for(int repeat = 0; repeat < 2000; repeat++) //use 2k here for NIF poster
                    {
                        if(uni() < (ProbAbsorb)*0.019*0.5) //uni() < (ProbAbsorb)*0.019*0.5
                        {
                            float AbsorptionLength = (-1.0f / AbsorbCoeff) * log(uni()); //how far xray travelled into material
                            while(AbsorptionLength > PathLength)
                            {
                                //need to guarantee that the photon is absorbed within the crystal.
                                //This is a stupid way to do this... If the absorption probablility is very low there's a very real chance of hanging here for a very long time.
                                AbsorptionLength = (-1.0f / AbsorbCoeff) * log(uni());
                            }
                            //isn't this caculation equivalent to checking if the photon is absorbed in the first place?
                            //should check if it is quicker just to do this test.

                            float ProbTransmit = exp( -1.0f * AbsorbCoeff * AbsorptionLength);
                            if(uni() < ProbTransmit)
                            {
                                //this checks if the photon would actually be transmitted (pass back out through the material)
                                continue;
                            }

                            double CosTheta = uni();//uni()*2.0 - 1.0;
                            if(CosTheta < DirectionMinCosTheta || CosTheta > DirectionMaxCosTheta)
                            {
                                continue;
                            }

                            double Phi = uni()*2.0*PI - PI; // uni()*2.0*PI;

                            if(Phi < DirectionMinPhi || Phi > DirectionMaxPhi)
                            {
                                continue;
                            }

                            Vector EmitDirection( acos(CosTheta), Phi, true, true);

                            /*if(EmitDirection.z <= 0.0)
                            {
                                cout << "Skipped" << endl;
                                continue;
                            }*/

                            /*if( EmitDirection.y > DirectionMaxY || EmitDirection.y < DirectionMinY ||
                                    EmitDirection.x > DirectionMaxX || EmitDirection.x < DirectionMinX)
                            {
                                continue;
                            }*/

                            Vector NewSource(Source.x - (Source.z/Direction.z)*Direction.x,
                                             Source.y - (Source.z/Direction.z)*Direction.y,
                                             0.0f);

                            //EmitDirection = EmitDirection.Normalized();

                            double CCDIntersectX = NewSource.x + ((CCDZ-NewSource.z)/(EmitDirection.z))*EmitDirection.x;
                            double CCDIntersectY = NewSource.y + ((CCDZ-NewSource.z)/(EmitDirection.z))*EmitDirection.y;

                            if( CCDIntersectX > 110.0 && CCDIntersectX < 115.0) //120-120.1
                            {
                                if( CCDIntersectY < 1.0 && CCDIntersectY > -1.0)
                                {
                                    FluoResults << CCDIntersectX << "\t" << CCDIntersectY << endl;
                                    AdvFluoResults << NewSource.x << "\t" << NewSource.y << "\t" <<
                                    EmitDirection.x << "\t" << EmitDirection.y << "\t" << EmitDirection.z << "\t" << endl;
                                    nFluoresced++;
                                }
                            }
                        }

                        else if(uni() < ProbScatter) //R < ProbScatter
                        {
                            float RockingCurve;
                            float BraggAngle = PD.PickBraggScatteringAngle( Energy, RockingCurve);

                            if(BraggAngle < 0.0f)
                            {
                                continue;
                            }

                            float PhiAngle = (uni() * 2.0*PI);//(UniformRand() * 1.0f) + PI-0.5f;

                            float ScatterAngle = (normal()*RockingCurve + BraggAngle)*2.0; //shift from standard normal to bragg peak parameters.

                            //Vector ScatterDirection = ScatterUnitVector( Direction, BraggAngle*2.0, PhiAngle); //Use this for no rocking curve.
                            Vector ScatterDirection = ScatterUnitVector( Direction, ScatterAngle, PhiAngle);

                            if(ScatterDirection.z <= 0.0f)
                            {
                                continue;
                            }

                            if( ScatterDirection.y > DirectionMaxY || ScatterDirection.y < DirectionMinY ||
                                    ScatterDirection.x > DirectionMaxX || ScatterDirection.x < DirectionMinX)
                            {
                                continue;
                            }

                            Vector NewSource(Source.x - (Source.z/Direction.z)*Direction.x,
                                             Source.y - (Source.z/Direction.z)*Direction.y,
                                             0.0f);

                            double CCDIntersectX = NewSource.x + ((CCDZ-NewSource.z)/(ScatterDirection.z))*ScatterDirection.x; //60
                            double CCDIntersectY = NewSource.y + ((CCDZ-NewSource.z)/(ScatterDirection.z))*ScatterDirection.y; //60

                            if( CCDIntersectX > 110.0f && CCDIntersectX < 115.0f)
                            {
                                if( CCDIntersectY < 1.0f && CCDIntersectY > -1.0f)
                                {
                                    DiffractResults << CCDIntersectX << "\t" << CCDIntersectY << "\t" << Energy << "\t" << BraggAngle << endl;
                                    AdvDiffractResults << NewSource.x << "\t" << NewSource.y << "\t" <<
                                    ScatterDirection.x << "\t" << ScatterDirection.y << "\t" << ScatterDirection.z << "\t" <<
                                    Energy << endl;
                                    nDiffracted++;
                                }
                            }
                        }
#ifdef FLU_ENABLE
                        //using R < blahblah gives clustering??
                        //Tantalum Fluo: http://www.nist.gov/data/PDFfiles/jpcrd473.pdf
                        //0.5 is from only creating x-rays that point upwards.
                        //else
#endif
                    }
                }
            }

            float Progress = float(EnergyTick)/float(nEPoints);

            if(Progress > ProgressCounter)
            {
                int ElapsedTime = time(NULL) - StartTime;
                cout << "Progress: " << Progress*100.0f << "%\t E = " << Energy << "\t(" << ElapsedTime << "s)"<< endl;
                ProgressCounter += 0.01f;
            }
        }
    }

    DiffractResults.close();
    FluoResults.close();

    AdvDiffractResults.close();
    AdvFluoResults.close();

    return 0;
}

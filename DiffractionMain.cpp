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
#include "SpectrumData.h"
#include "MPIUtils.h"

#include "mpi.h"

using namespace std;

#ifndef base_generator_type
//typedef boost::mt11213b base_generator_type;
typedef boost::mt19937 base_generator_type; //Docs say that this is marginally slower but I'm unable to measure a repeatable difference. Supposedly this also has better random properties?
#endif



struct EnergyWorkUnit
{
    unsigned long int StartEnergyTick, EndEnergyTick;
};


//http://www.boost.org/doc/libs/1_46_1/libs/random/example/random_demo.cpp
//http://www.boost.org/doc/libs/1_52_0/doc/html/boost_random/reference.html


#define FLUO_DISABLE
//#define FORCE_DIFFRACTION

int main(int argc, char *argv[])
{
    
    
    
    MPI_Init(&argc,&argv);
        
    int NumProcessors, ProcessorId;
    MPI_Comm_size(MPI_COMM_WORLD,&NumProcessors);
    MPI_Comm_rank(MPI_COMM_WORLD,&ProcessorId);
    
#ifdef FLUO_DISABLE
    
    cout << "Fluorscence disabled! Diffraction forced on" << endl;;
    
#endif
    
    ifstream datafile("InputScript.txt");
    if(datafile.is_open() == false)
    {
        cout << "Error: Failed to open InputScript.txt" << endl;
        exit(1);
    }
    std::map<std::string,std::string> InputData;
    AddToMapFromFile(datafile, InputData);
    datafile.close();

    
    CCD CCDCamera = GenerateCCDFromInputScript("InputScript.txt");

    Vector CCDCorners[4];

    CCDCamera.GenerateCCDCorners( CCDCorners[0], CCDCorners[1], CCDCorners[2], CCDCorners[3]);

    double CCDXMin = CCDCorners[0].x, CCDXMax = CCDCorners[0].x, CCDYMin = CCDCorners[0].y, CCDYMax = CCDCorners[0].y;
    
    cout << "CCDBounds:" << endl;
    
    
    for(int i = 0; i<4; i++)
    {
        cout << i << ":\t"; CCDCorners[i].Print();
        if(CCDCorners[i].x < CCDXMin)
        {
            CCDXMin = CCDCorners[i].x;
        }
        if(CCDCorners[i].y < CCDYMin)
        {
            CCDYMin = CCDCorners[i].y;
        }
        if(CCDCorners[i].x > CCDXMax)
        {
            CCDXMax = CCDCorners[i].x;
        }
        if(CCDCorners[i].y > CCDYMax)
        {
            CCDYMax = CCDCorners[i].y;
        }

    }
   
    Vector CrystalOrigin(0,-0.1,0);
    double CrystalXLength = 0.1;
    double CrystalYLength = 0.2;
    
    VectorFromMap("CrystalOrigin",InputData, CrystalOrigin);
    DoubleFromMap("CrystalXLength", InputData, CrystalXLength);
    DoubleFromMap("CrystalYLength", InputData, CrystalYLength);
    
    cout << "CrystalOrigin:\t"; CrystalOrigin.Print();
    cout << "CrystalDimensions:\tX:\t" << CrystalXLength << "\tY:\t" << CrystalYLength << endl;
    
    Vector CrystalCorners[4];
    
    CrystalCorners[0] = Vector(CrystalOrigin);
    CrystalCorners[1] = Vector(CrystalOrigin.x,CrystalOrigin.y+CrystalYLength,0);
    CrystalCorners[2] = Vector(CrystalOrigin.x+CrystalXLength,CrystalOrigin.y,0);
    CrystalCorners[3] = Vector(CrystalOrigin.x+CrystalXLength,CrystalOrigin.y+CrystalYLength,0);
    
    cout << "Crystal Corners:" << endl;
    for(int i = 0; i<4; i++)
    {
        CrystalCorners[i].Print();
    }
    
    
    double DirectionMinCosTheta = 0.0, DirectionMaxCosTheta = 0.0, DirectionMinPhi = 0.0, DirectionMaxPhi = 0.0;
    double DirectionMinX = 5, DirectionMinY = 5, DirectionMaxX = -5, DirectionMaxY = -5;


    bool bFirstIteration = true;

    //loop through each corner of the crystal, then loop through each of the CCD corners to find the min/max for the solid angle rejection parameters.
    for(int iCrystalCorner = 0; iCrystalCorner<4; iCrystalCorner++)
    {
        for( int iCCDCorner = 0; iCCDCorner<4; iCCDCorner++)
        {
            Vector Direction = CCDCorners[iCCDCorner] - CrystalCorners[iCrystalCorner];
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


    double SolidAngleTolerance = 0.005; //% tolerance for solid angle rejection. (how much we add/take away to the x&y limits)

    DoubleFromMap("SolidAngleTolerance", InputData, SolidAngleTolerance);

    DirectionMinX -= DirectionMinX*SolidAngleTolerance;
    DirectionMinY -= DirectionMinY*SolidAngleTolerance;
    DirectionMaxX += DirectionMaxX*SolidAngleTolerance;
    DirectionMaxY += DirectionMaxY*SolidAngleTolerance;

    DirectionMinCosTheta -= DirectionMinCosTheta*SolidAngleTolerance;
    DirectionMinPhi -= DirectionMinPhi*SolidAngleTolerance;
    DirectionMaxCosTheta += DirectionMaxCosTheta*SolidAngleTolerance;
    DirectionMaxPhi += DirectionMaxPhi*SolidAngleTolerance;

    cout << "Emitted ray bounds:" << endl;
    cout << "X:\t" << DirectionMinX << "\t" << DirectionMaxX << endl;
    cout << "Y:\t" << DirectionMinY << "\t" << DirectionMaxY << endl;
    
    
    cout << "CosTheta:\t" << DirectionMinCosTheta << "\t" << DirectionMaxCosTheta << endl;
    cout << "Phi:\t" << DirectionMinPhi << "\t" << DirectionMaxPhi << endl;


    base_generator_type generator(42580*(ProcessorId+500));
    
    
    boost::uniform_real<> uni_dist(0,1);
    boost::normal_distribution<> normal_dist(0,1);
    boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
    boost::variate_generator<base_generator_type&, boost::normal_distribution<> > normal(generator, normal_dist);
    

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
    //ofstream DiffractResults( "DiffractResults.txt" );
    //ofstream FluoResults( "FluoResults.txt" );
    ofstream DiffractResults( CreateProcessorUniqueFilename("DiffractResults", ".txt") );
    ofstream FluoResults( CreateProcessorUniqueFilename("FluoResults", ".txt") );
    

    //full ray information
    //ofstream AdvDiffractResults( "AdvDiffractResults.txt");
    //ofstream AdvFluoResults( "AdvFluoResults.txt" );
    ofstream AdvDiffractResults( CreateProcessorUniqueFilename("AdvDiffractResults", ".txt"));
    ofstream AdvFluoResults( CreateProcessorUniqueFilename("AdvFluoResults", ".txt") );

    double CCDZ = CCDCamera.CCDOrigin.z;

    //DoubleFromMap("SourceZ", InputData, SourceZ);
    //DoubleFromMap("CCDZ", InputData, CCDZ);

    Vector Source( -3.0f, 0.0f, 5.0);    
    VectorFromMap("Source", InputData, Source);
    double SourceZ = Source.z;

    double MinE = 3.0;//4.23;
    double MaxE = 9.0;//4.28;
    double DeltaE = 0.0005;//0.0005f;

    DoubleFromMap("MinEnergy", InputData, MinE);
    DoubleFromMap("MaxEnergy", InputData, MaxE);
    DoubleFromMap("DeltaE", InputData, DeltaE);

    int nEPoints = (MaxE - MinE)/DeltaE;
           
    
    double minSourceDirectionX = (CrystalOrigin - Source).x;
    double maxSourceDirectionX = minSourceDirectionX + CrystalXLength;
    double DeltaSourceDirectionX        = 0.01;
    DoubleFromMap("TempSourceDeltaX", InputData, DeltaSourceDirectionX);
    
    if(minSourceDirectionX > maxSourceDirectionX)
    {
        swap(minSourceDirectionX,maxSourceDirectionX);
    }
    
    cout << "X Source Ranges:\t" << minSourceDirectionX << "\t" << maxSourceDirectionX << endl;
    
    
    double minSourceDirectionY = (CrystalOrigin - Source).y;
    double maxSourceDirectionY = minSourceDirectionY + CrystalYLength;
    double DeltaSourceDirectionY        = 0.01;
    DoubleFromMap("TempSourceDeltaY", InputData, DeltaSourceDirectionY);

    
    if(minSourceDirectionY > maxSourceDirectionY)
    {
        swap(minSourceDirectionY,maxSourceDirectionY);
    }
    
    cout << "Y Source Ranges:\t" << minSourceDirectionY << "\t" << maxSourceDirectionY << endl;
    
    
    int nDirectionXPoints = (maxSourceDirectionX - minSourceDirectionX)/DeltaSourceDirectionX; 
    int nDirectionYPoints = (maxSourceDirectionY - minSourceDirectionY)/DeltaSourceDirectionY;
       
    
    int nRepeats = 500;
    IntFromMap("Repeats", InputData, nRepeats);
    cout << nRepeats << " Repeats" << endl;
    
    long unsigned int nPhotons = 0;
    
    bool bRockingCurve = true;
    int iRockingCurve;
    IntFromMap("RockingCurve",InputData,iRockingCurve);
    
    bRockingCurve = int(iRockingCurve);
    if(!bRockingCurve)
    {
        cout << "Rocking Curve Disabled" << endl;
    }
    
    
    std::string SpectrumDataFilename;
    
    StringFromMap("SpectrumData", InputData, SpectrumDataFilename);
    SpectrumData Spectrum( MinE, MaxE, 5000 );
    Spectrum.LoadData(SpectrumDataFilename.c_str());
    
    cout << "Generating Energy Intervals" << endl;
    std::vector<std::pair< double, double >> Intervals;
    
    if(ProcessorId == 0)
    {    
        Spectrum.GetEnergyIntervals( NumProcessors, &Intervals);
        if(Intervals.size() != NumProcessors)
        {
            cout << "Spectrum Integration failed! Quitting!" << endl;
            exit(1);
        }
    }

    
    int nDiffracted = 0, nFluoresced = 0;
    
    long unsigned int StartTime = time(NULL);
    

    if(ProcessorId == 0)
    {
        unsigned long int WorkSize = nEPoints/NumProcessors;
        
        for(unsigned long int i=0;i<NumProcessors;i++)
        {
            EnergyWorkUnit Work;
            //Work.StartEnergyTick = i * WorkSize;
            //Work.EndEnergyTick = (i+1) * WorkSize;
            
            Work.StartEnergyTick = (Intervals[i].first - MinE)/DeltaE;
            Work.EndEnergyTick = (Intervals[i].second - MinE)/DeltaE;
            
            
            if( i == (NumProcessors - 1) ) //do this so we get the final energy loop for the last core
            {
                Work.EndEnergyTick++;
            }
            MPI_Send(&Work, sizeof(EnergyWorkUnit), MPI_BYTE, int(i), 0, MPI_COMM_WORLD);
        }
    }
    
    MPI_Status Stat;
    EnergyWorkUnit WorkToDo;    
    MPI_Recv(&WorkToDo, sizeof(EnergyWorkUnit), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &Stat);
    cout << ProcessorId << ":\t" << WorkToDo.StartEnergyTick << "  -->  " << WorkToDo.EndEnergyTick << endl;
    
    
    float ProgressCounter = 0.0f;    
    float AbsorbCoeffFluo = TaMuData.GetAbsorbCoeffDataPoint(EnergyToWavelength(1.710));
    
    
    float MinTheta = Deg2Rad( 0 );
    float MaxTheta  = Deg2Rad( +2.862 );
    float MinPhi = Deg2Rad( 0 );
    float MaxPhi  = Deg2Rad( 360 );
    
    float MinCosTheta = cos( MinTheta );
    float MaxCosTheta = cos( MaxTheta );
    
    if (MinCosTheta > MaxCosTheta)
    {
        swap( MinCosTheta, MaxCosTheta);
    }
    
    if (MinPhi > MaxPhi)
    {
        swap( MinPhi, MaxPhi);
    }
    
    float deltaCosTheta = (MaxCosTheta-MinCosTheta)/200.0;
    float deltaPhi = (MaxPhi-MinPhi)/100.0;
    
    Vector SourceToCrystal( 3.830, 0, -3.214);    
    //Direction = TransformToNewFrame(Direction, SourceToCrystal.GetTheta(), SourceToCrystal.GetPhi());
    
    for( unsigned long int EnergyTick = WorkToDo.StartEnergyTick; EnergyTick < WorkToDo.EndEnergyTick; EnergyTick++)
    {
        double Energy = MinE + DeltaE*EnergyTick;
        
        double RelativeSourceIntensity = Spectrum.GetSpectrumDataPoint(Energy);
        int CorrectedRepeats = floor(RelativeSourceIntensity*double(nRepeats) + 0.5); //round to nearest integer
        
        float ProbScatter = PD.GetModifiedScatterProb(Energy);
        float AbsorbCoeff = TaMuData.GetAbsorbCoeffDataPoint( EnergyToWavelength(Energy) );
          
        
        //for(int iXDirectionCounter = 0; iXDirectionCounter <= nDirectionXPoints; iXDirectionCounter++)
        for( float CosTheta = MinCosTheta; CosTheta <= MaxCosTheta; CosTheta += deltaCosTheta)
        {
            //double x = minSourceDirectionX + DeltaSourceDirectionX*iXDirectionCounter;
            
            //for(int iYDirectionCounter = 0; iYDirectionCounter <= nDirectionYPoints; iYDirectionCounter++)
            for( float Phi = MinPhi; Phi <= MaxPhi; Phi += deltaPhi)
            {
                //double y = minSourceDirectionY + DeltaSourceDirectionY*iYDirectionCounter;
                //construct a vector that goes down SourceZ units and across a variable amount
                //Vector Direction( x, y, -1.0*SourceZ);
                //Direction = Direction.Normalized();
                
                
                Vector Direction( acos(CosTheta), Phi, true, true);
                Direction = TransformToNewFrame(Direction, SourceToCrystal.GetTheta(), SourceToCrystal.GetPhi());
                
                float PathLength = fabs(Thickness/Direction.z); //length xray takes through crystal
                double ProbAbsorb = 1.0 - exp( -1.0 * AbsorbCoeff * PathLength);
                
                //float CorrectedProbScatter = ProbScatter;
                
                /*if( (1.0 - ProbAbsorb) > 0.0)
                 {
                 CorrectedProbScatter = ProbScatter/(1.0-ProbAbsorb);
                 } */
                
                ProbAbsorb = ProbAbsorb/(1.0-ProbScatter); //Is this negligble? Is this even correct??
                
                
                for(int repeat = 0; repeat < CorrectedRepeats; repeat++) //use 2k here for NIF poster
                {
                    nPhotons++;
#ifndef FORCE_DIFFRACTION
                    if(uni() < ProbScatter)
#else
                    if(1)
#endif
                    {
                        float RockingCurve;
                        float BraggAngle = PD.PickBraggScatteringAngle( Energy, RockingCurve);
                        
                        if(BraggAngle < 0.0f)
                        {
                            continue;
                        }
                        
                        float PhiAngle = (uni() * 2.0*PI);//(UniformRand() * 1.0f) + PI-0.5f;
                        
                        float ScatterAngle;
                        
                        
                        if(bRockingCurve)
                        {
                            //width is for peak at 2theta (not at bragg angle) so no need to multiply it by 2.
                            ScatterAngle = normal()*RockingCurve + BraggAngle*2.0; //shift from standard normal to bragg peak parameters.
                        }
                        else
                        {
                            ScatterAngle = BraggAngle*2.0;
                        }
                        
                        Vector ScatterDirection = ScatterUnitVector( Direction, ScatterAngle, PhiAngle);
                        
                        if(ScatterDirection.z <= 0.0f)
                        {
                            continue;
                        }
                        
                        float CosTheta = cos(ScatterDirection.GetTheta());
                        float Phi = ScatterDirection.GetPhi();
                        
                        /*if( ScatterDirection.y > DirectionMaxY || ScatterDirection.y < DirectionMinY ||
                         ScatterDirection.x > DirectionMaxX || ScatterDirection.x < DirectionMinX)
                         {
                         continue;
                         }*/
                        
                        if(CosTheta < DirectionMinCosTheta || CosTheta > DirectionMaxCosTheta)
                        {
                            continue;
                        }
                        
                        if(Phi < DirectionMinPhi || Phi > DirectionMaxPhi)
                        {
                            continue;
                        }
                        
                        Vector NewSource(Source.x - (Source.z/Direction.z)*Direction.x,
                                         Source.y - (Source.z/Direction.z)*Direction.y,
                                         0.0f);
                        
                        double CCDIntersectX = NewSource.x + ((CCDZ-NewSource.z)/(ScatterDirection.z))*ScatterDirection.x; //60
                        double CCDIntersectY = NewSource.y + ((CCDZ-NewSource.z)/(ScatterDirection.z))*ScatterDirection.y; //60
                        
                        /*if( CCDIntersectX > CCDXMin && CCDIntersectX < CCDXMax)
                         {
                         if( CCDIntersectY < CCDYMax && CCDIntersectY > CCDYMin)
                         {
                         DiffractResults << CCDIntersectX << "\t" << CCDIntersectY << "\t" << Energy << "\t" << BraggAngle << endl;
                         AdvDiffractResults << NewSource.x << "\t" << NewSource.y << "\t" <<
                         ScatterDirection.x << "\t" << ScatterDirection.y << "\t" << ScatterDirection.z << "\t" <<
                         Energy << endl;
                         nDiffracted++;
                         }
                         }*/
                        DiffractResults << CCDIntersectX << "\t" << CCDIntersectY << "\t" << Energy << "\t" << BraggAngle << endl;
                        AdvDiffractResults << NewSource.x << "\t" << NewSource.y << "\t" <<
                        ScatterDirection.x << "\t" << ScatterDirection.y << "\t" << ScatterDirection.z << "\t" <<
                        Energy << endl;
                        nDiffracted++;
                    }
#ifndef FLUO_DISABLE
                    //Tantalum Fluo: http://www.nist.gov/data/PDFfiles/jpcrd473.pdf
                    //0.5 is from only creating x-rays that point upwards.
                    else if(uni() < (ProbAbsorb)*0.019*0.5)//if(uni() < (ProbAbsorb))//
                    {
                        /*if(uni() > 0.019*0.5)
                        {
                            continue;
                        }*/
                        float AbsorptionLength = (-1.0f / AbsorbCoeff) * log(uni()); //how far xray travelled into material
                        /*while(AbsorptionLength > PathLength)
                        {
                            //need to guarantee that the photon is absorbed within the crystal.
                            //This is a stupid way to do this... If the absorption probablility is very low there's a very real chance of hanging here for a very long time.
                            AbsorptionLength = (-1.0f / AbsorbCoeff) * log(uni());
                        }*/
                        //isn't this caculation equivalent to checking if the photon is absorbed in the first place?
                        //should check if it is quicker just to do this test.
                                                
                        if(AbsorptionLength > PathLength)
                        {
                            //this almost never gets hit for reasonably sized crystals                            
                            AbsorptionLength = PathLength;
                        }
                        
                        float ProbTransmit = exp( -1.0f * AbsorbCoeffFluo * AbsorptionLength);
                        
                        if(uni() > ProbTransmit) //if absorbed, continue to next photon
                        {
                            //this checks if the photon would actually be transmitted (pass back out through the material)
                            continue;
                        }
                        
                        double CosTheta = uni();//uni()*2.0 - 1.0;
                        if(CosTheta < DirectionMinCosTheta || CosTheta > DirectionMaxCosTheta)
                        {
                            continue;
                        }
                        
                        //need to change phi interval to [pi,-pi] so that it passes the phi limits below
                        //which are created in that same interval
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
                        
                        /*if( CCDIntersectX > CCDXMin && CCDIntersectX < CCDXMax)
                         {
                         if( CCDIntersectY < CCDYMax && CCDIntersectY > CCDYMin)
                         {
                         FluoResults << CCDIntersectX << "\t" << CCDIntersectY << endl;
                         AdvFluoResults << NewSource.x << "\t" << NewSource.y << "\t" <<
                         EmitDirection.x << "\t" << EmitDirection.y << "\t" << EmitDirection.z << "\t" << endl;
                         nFluoresced++;
                         }
                         }*/
                        FluoResults << CCDIntersectX << "\t" << CCDIntersectY << endl;
                        AdvFluoResults << NewSource.x << "\t" << NewSource.y << "\t" <<
                        EmitDirection.x << "\t" << EmitDirection.y << "\t" << EmitDirection.z << "\t" << endl;
                        nFluoresced++;
                    }
#endif
                }
            }
        }
        
        float Progress = float(EnergyTick - WorkToDo.StartEnergyTick)/float(WorkToDo.EndEnergyTick - WorkToDo.StartEnergyTick);
        
        if(Progress > ProgressCounter)
        {
            int ElapsedTime = int(time(NULL) - StartTime);
            printf( "Proc %d: Progress: %.0f%%\tE = %.3f\t(%is)\n", ProcessorId, Progress*100.0, Energy, ElapsedTime); //output looks nicer this way
            //cout << "Proc " << ProcessorId << ": Progress: " << Progress*100.0f << "%\t E = " << Energy << "\t(" << ElapsedTime << "s)"<< endl;
            ProgressCounter += 0.01f;
        }
    }
    

    DiffractResults.close();
    FluoResults.close();

    AdvDiffractResults.close();
    AdvFluoResults.close();
    
    
    cout << "Proc " << ProcessorId << ": " <<  nPhotons    << " Photons"   << endl;
    cout << "Proc " << ProcessorId << ": " <<  nDiffracted << " Diffracted photons" << endl;
    cout << "Proc " << ProcessorId << ": " <<  nFluoresced << " Fluoresced photons" << endl;
    
    cout << "Done!" << endl;

    MPI_Barrier(MPI_COMM_WORLD);
    
    //cat DiffractResults_P0.txt DiffractResults_P1.txt DiffractResults_P2.txt DiffractResults_P3.txt > DiffractResults.txt
    if( ProcessorId == 0)
    {
        system(CreateConcatCommand("DiffractResults", ".txt").c_str());
        system(CreateConcatCommand("AdvDiffractResults", ".txt").c_str());
        system(CreateConcatCommand("FluoResults", ".txt").c_str());
        system(CreateConcatCommand("AdvFluoResults", ".txt").c_str());
    }

    MPI_Finalize();
    
    
    
    return 0;
}

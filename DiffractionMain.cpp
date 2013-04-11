#include "mpi.h"
#include "MPIUtils.h"

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
#include "FileReading.h"
#include "SpectrumData.h"
#include "FluorescenceData.h"

#include "CCD.h"

using namespace std;

#ifndef base_generator_type
//typedef boost::mt11213b base_generator_type;
typedef boost::mt19937 base_generator_type; //Docs say that this is marginally slower but I'm unable to measure a repeatable difference. Supposedly this also has better random properties?
#endif



struct EnergyWorkUnit
{
    unsigned int StartEnergyTick, EndEnergyTick;
};





//#define FLUO_DISABLE
//#define FORCE_DIFFRACTION

int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);
    MPI_cout << endl;
    MPI_cout << "////////////////////////////////////////" << endl;
    MPI_cout << "//Now running Main Diffraction Program//" << endl;
    MPI_cout << "////////////////////////////////////////" << endl;
    MPI_cout << endl;
        
    int NumProcessors, ProcessorId;
    MPI_Comm_size(MPI_COMM_WORLD,&NumProcessors);
    MPI_Comm_rank(MPI_COMM_WORLD,&ProcessorId);
    
	
	MPI_cout << "Num Processors: " << NumProcessors << endl;
	//cout << "Processor " <<	ProcessorId << " online" << endl;
	
#ifdef FLUO_DISABLE    
    MPI_cout << "Fluorscence disabled!" << endl;    
#endif
    
#ifdef FORCE_DIFFRACTION    
    MPI_cout << "Diffraction forced on!" << endl;    
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
    
    MPI_cout << "CCDBounds:" << endl;
        
    for(int i = 0; i<4; i++)
    {
        MPI_cout << i << ":\t"; if(ProcessorId == 0){CCDCorners[i].Print();}
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
    
    MPI_cout << "CrystalOrigin:\t";if(ProcessorId == 0){CrystalOrigin.Print();}
    MPI_cout << "CrystalDimensions:\tX:\t" << CrystalXLength << "\tY:\t" << CrystalYLength << endl;
    

    Vector Source( -3.0f, 0.0f, 5.0);
    VectorFromMap("Source", InputData, Source);
    
    double DirectionMinCosTheta = 0.0, DirectionMaxCosTheta = 0.0, DirectionMinPhi = 0.0, DirectionMaxPhi = 0.0;
     
    Vector SourceToOrigin = -1.0*Source;
    
    double SourceDivergence = 0.0;
    DoubleFromMap("SourceDivergence", InputData, SourceDivergence);
    
    
    float MinTheta = 0.0;
    float MaxTheta  = Deg2Rad( +SourceDivergence );
    float MinPhi = 0.0;
    float MaxPhi  = 2.0*PI;
    
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
        
    int NumThetaSteps = 1;
    int NumPhiSteps = 1;
    IntFromMap( "NumThetaSteps", InputData, NumThetaSteps);
    IntFromMap( "NumPhiSteps", InputData, NumPhiSteps);
    
    double deltaCosTheta = (MaxCosTheta-MinCosTheta)/double(NumThetaSteps);
    double deltaPhi = (MaxPhi-MinPhi)/double(NumPhiSteps);
    
    Vector SourceToCrystal = -1.0*Source;
    double SourceToCrystalTheta = SourceToCrystal.GetTheta();
    double SourceToCrystalPhi = SourceToCrystal.GetPhi();
    
    bool bFirstIteration = true;

    for( int CosThetaStep = 0; CosThetaStep <= NumThetaSteps; CosThetaStep ++)
    {
        double SourceCosTheta = MinCosTheta + double(CosThetaStep)*deltaCosTheta;
        
        for( int PhiStep = 0; PhiStep <= NumPhiSteps; PhiStep ++)
        {
            double SourcePhi = MinPhi + double(PhiStep)*deltaPhi;
            
            Vector SourceToCrystal( acos(SourceCosTheta), SourcePhi, true, true);
            SourceToCrystal = TransformToNewFrame(SourceToCrystal, SourceToOrigin.GetTheta(), SourceToOrigin.GetPhi());
            SourceToCrystal = SourceToCrystal.Normalized();
            
            Vector CrystalIntersection(Source.x - (Source.z/SourceToCrystal.z)*SourceToCrystal.x,
                                       Source.y - (Source.z/SourceToCrystal.z)*SourceToCrystal.y,
                                       0.0f); //intersection of ray with z=0 plane
                        
            for( int iCCDCorner = 0; iCCDCorner<4; iCCDCorner++)
            {
                Vector Direction = CCDCorners[iCCDCorner] - CrystalIntersection;
                Direction = Direction.Normalized();
                
                double CosTheta = cos(Direction.GetTheta());
                double Phi = Direction.GetPhi();
                
                if(bFirstIteration)
                {                    
                    DirectionMinCosTheta = CosTheta;
                    DirectionMaxCosTheta = CosTheta;
                    DirectionMinPhi = Phi;
                    DirectionMaxPhi = Phi;
                    
                    bFirstIteration = false;
                    continue;
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
    }


    double SolidAngleTolerance = 0.005; //% tolerance for solid angle rejection. (how much we add/take away to the x&y limits)

    DoubleFromMap("SolidAngleTolerance", InputData, SolidAngleTolerance);

    DirectionMinCosTheta -= DirectionMinCosTheta*SolidAngleTolerance;
    DirectionMinPhi -= DirectionMinPhi*SolidAngleTolerance;
    DirectionMaxCosTheta += DirectionMaxCosTheta*SolidAngleTolerance;
    DirectionMaxPhi += DirectionMaxPhi*SolidAngleTolerance;

    MPI_cout << "Emitted ray bounds:" << endl;
    
    
    MPI_cout << "CosTheta:\t" << DirectionMinCosTheta << "\t" << DirectionMaxCosTheta << endl;
    MPI_cout << "Phi:\t" << DirectionMinPhi << "\t" << DirectionMaxPhi << endl;

	
	
	//http://www.boost.org/doc/libs/1_46_1/libs/random/example/random_demo.cpp
	//http://www.boost.org/doc/libs/1_52_0/doc/html/boost_random/reference.html
	
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

    
    
    std::string MuDataFilename;
    StringFromMap("AbsorptionData", InputData, MuDataFilename);
        
    double MinE = 3.0;//4.23;
    double MaxE = 9.0;//4.28;
    double DeltaE = 0.0005;//0.0005f;
    
    DoubleFromMap("MinEnergy", InputData, MinE);
    DoubleFromMap("MaxEnergy", InputData, MaxE);
    DoubleFromMap("DeltaE", InputData, DeltaE);
    	
	
    double MinWavelength = EnergyToWavelength(MaxE); MinWavelength -= MinWavelength*0.02;
    double MaxWavelength = EnergyToWavelength(MinE);
    
    if( MaxWavelength < EnergyToWavelength(0.4) )
    {
        MaxWavelength = EnergyToWavelength(0.4);
    }
    
    MaxWavelength += MaxWavelength*0.02;
    
    AbsorbCoeffData MuData( MinWavelength, MaxWavelength, 1, 10000, MuDataFilename.c_str());


#ifndef FLUO_DISABLE
	std::string EmissionLineData;
	std::string ShellProbData;
	
	StringFromMap("EmissionLineData", InputData, EmissionLineData);
	StringFromMap("ShellProbabilityData", InputData, ShellProbData);
		
	FluorescenceData FluoData( MinE, MaxE, 2000, EmissionLineData.c_str(), ShellProbData.c_str(), &uni );
#endif
	


	
    double Thickness = 500000; //50 microns in A
    DoubleFromMap("CrystalThickness", InputData, Thickness);

	
    //simple detector at X/Y plane
    ofstream DiffractResults( CreateProcessorUniqueFilename("DiffractResults", ".txt").c_str() );
    ofstream FluoResults( CreateProcessorUniqueFilename("FluoResults", ".txt").c_str() );

    //full ray information
    ofstream AdvDiffractResults( CreateProcessorUniqueFilename("AdvDiffractResults", ".txt").c_str() );
    ofstream AdvFluoResults( CreateProcessorUniqueFilename("AdvFluoResults", ".txt").c_str() );

    double CCDZ = CCDCamera.CCDOrigin.z;


    int nRepeats = 500;
    IntFromMap("Repeats", InputData, nRepeats);
    MPI_cout << nRepeats << " Repeats" << endl;
    
    long unsigned int nPhotons = 0;
    
    bool bRockingCurve = true;
    int iRockingCurve;
    IntFromMap("RockingCurve",InputData,iRockingCurve);
    
    bRockingCurve = int(iRockingCurve);
    if(!bRockingCurve)
    {
        MPI_cout << "Rocking Curve Disabled" << endl;
    }
    

	
    std::string SpectrumDataFilename;    
    StringFromMap("SpectrumData", InputData, SpectrumDataFilename);
    SpectrumData Spectrum( MinE, MaxE, 1, 5000, SpectrumDataFilename.c_str() );
        
    int nDiffracted = 0, nFluoresced = 0;
    
    MPI_cout << "Generating Energy Intervals" << endl;
	
	EnergyWorkUnit WorkToDo;
	
    if(ProcessorId == 0)
    {
        
        std::vector<std::pair< double, double > > Intervals;
        Spectrum.GetEnergyIntervals( NumProcessors, &Intervals);
		cout << "Finished generating intervals!" << endl;
        if(Intervals.size() != NumProcessors)
        {
            cout << "Spectrum Integration failed! Quitting!" << endl;
            exit(1);
        }
        cout << "Sending intervals...";
		
		EnergyWorkUnit WorkUnits[NumProcessors];
		
        for(int i = 0;i<NumProcessors;i++)
        {
            //EnergyWorkUnit Work;
            
            WorkUnits[i].StartEnergyTick = (Intervals[i].first - MinE)/DeltaE;
            WorkUnits[i].EndEnergyTick = (Intervals[i].second - MinE)/DeltaE;
                        
            if( i == (NumProcessors - 1) ) //do this so we get the final energy loop for the last core
            {
                WorkUnits[i].EndEnergyTick++;
            }
			MPI_Request R;
			if( i != 0 )
			{
				MPI_Isend(&(WorkUnits[i]), sizeof(EnergyWorkUnit), MPI_BYTE, int(i), 0, MPI_COMM_WORLD, &R);
			}
			else
			{
				WorkToDo.StartEnergyTick = WorkUnits[0].StartEnergyTick;
				WorkToDo.EndEnergyTick = WorkUnits[0].EndEnergyTick;
			}
			//cout << "Sent Processor " << i << ": " << WorkUnits[i].StartEnergyTick << "-->" << WorkUnits[i].EndEnergyTick << endl;
        }
		cout << "Intervals sent" << endl;
    }
    
    MPI_Status Stat;
    
	MPI_cout << "Receiving intervals" << endl;
	
	
	for( int i = 0 ; i < NumProcessors; i++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		if( ProcessorId == i)
		{
			if( i != 0)
			{
				MPI_Recv(&WorkToDo, sizeof(EnergyWorkUnit), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &Stat);
			}
			cout << ProcessorId << ":\t" << WorkToDo.StartEnergyTick << "  -->  " << WorkToDo.EndEnergyTick << endl;
		}
    }
    
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	
    float ProgressCounter = 0.0f;
    


    
    long unsigned int StartTime = time(NULL);
    
    //double FluoYield = 0.019;
    
	//cout << "Processor " << ProcessorId << " starting" << endl;
	
    for( unsigned long int EnergyTick = WorkToDo.StartEnergyTick; EnergyTick < WorkToDo.EndEnergyTick; EnergyTick++)
    {
        double Energy = MinE + DeltaE*EnergyTick;
        
        double RelativeSourceIntensity = Spectrum.GetSpectrumDataPoint(Energy);
        double SpectrumCorrectedRepeats = RelativeSourceIntensity*double(nRepeats);//floor(RelativeSourceIntensity*double(nRepeats) + 0.5); //round to nearest integer
        
        float ProbScatter = PD.GetModifiedScatterProb(Energy);
        float AbsorbCoeff = MuData.GetAbsorbCoeffDataPoint( EnergyToWavelength(Energy) );
#ifndef FLUO_DISABLE
		FluoData.PreSelectEnergy(Energy); //pre-gather data as high up the loop as we can.
#endif
        for( int CosThetaStep = 0; CosThetaStep <= NumThetaSteps; CosThetaStep ++)
        {
            double CosTheta = MinCosTheta + double(CosThetaStep)*deltaCosTheta;
            
            double SolidAngleCorrectedRepeats = SpectrumCorrectedRepeats;
            if(CosThetaStep == 0 || CosThetaStep == NumThetaSteps)
            {
                SolidAngleCorrectedRepeats *= 0.5; //correction to account for first and last costheta steps hitting a smaller solid angle
            }
            
            int CorrectedRepeats = floor( SolidAngleCorrectedRepeats + 0.5 ); //floor could probably be replaced by an int cast since we're only dealing +ve numbers here.
            
            for( int PhiStep = 0; PhiStep < NumPhiSteps; PhiStep ++) //Do I want to go around by 2pi or slightly less than 2pi?
            {
                double Phi = MinPhi + double(PhiStep)*deltaPhi;
                
                Vector Direction( acos(CosTheta), Phi, true, true);
                Direction = TransformToNewFrame(Direction, SourceToCrystalTheta, SourceToCrystalPhi);                
                float PathLength = fabs(Thickness/Direction.z); //length xray takes through crystal
                
                double ProbAbsorb = 1.0 - exp( -1.0 * AbsorbCoeff * PathLength);
                ProbAbsorb *= 0.5; // accounts for xrays (that we care about) only being emitted in +ve Z
                ProbAbsorb += ProbScatter; //this correction is likely to be negligible
                
                //float CorrectedProbScatter = ProbScatter;
                
                /*if( (1.0 - ProbAbsorb) > 0.0)
                 {
                 CorrectedProbScatter = ProbScatter/(1.0-ProbAbsorb);
                 } */
                
                //ProbAbsorb = ProbAbsorb/(1.0-ProbScatter); //Is this negligble? Is this even correct??
                
                                
                for(int repeat = 0; repeat < CorrectedRepeats; repeat++)
                {
                    double R = uni();
                    nPhotons++;
#ifndef FORCE_DIFFRACTION
                    if(R < ProbScatter)
#else
                    if(true)
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
                        if(CosTheta < DirectionMinCosTheta || CosTheta > DirectionMaxCosTheta)
                        {
                            continue;
                        }
						
						
						float Phi = ScatterDirection.GetPhi();                        
                        if(Phi < DirectionMinPhi || Phi > DirectionMaxPhi)
                        {
                            continue;
                        }
                        
                        Vector NewSource(Source.x - (Source.z/Direction.z)*Direction.x,
                                         Source.y - (Source.z/Direction.z)*Direction.y,
                                         0.0f);
                        
                        double CCDIntersectX = NewSource.x + ((CCDZ-NewSource.z)/(ScatterDirection.z))*ScatterDirection.x; //60
                        double CCDIntersectY = NewSource.y + ((CCDZ-NewSource.z)/(ScatterDirection.z))*ScatterDirection.y; //60
                       

                        DiffractResults << CCDIntersectX << "\t" << CCDIntersectY << "\t" << Energy << "\t" << BraggAngle << endl;
                        AdvDiffractResults << NewSource.x << "\t" << NewSource.y << "\t" <<
                        ScatterDirection.x << "\t" << ScatterDirection.y << "\t" << ScatterDirection.z << "\t" <<
                        Energy << endl;
                        nDiffracted++;
                    }
#ifndef FLUO_DISABLE
                    //Tantalum Fluo: http://www.nist.gov/data/PDFfiles/jpcrd473.pdf
                    //0.5 is from only creating x-rays that point upwards.
                    else if(R < ProbAbsorb)
                    {
						//need to change phi interval to [pi,-pi] so that it passes the phi limits below
                        //which are created in that same interval
                        double Phi = uni()*2.0*PI - PI; // uni()*2.0*PI;
                        
                        if(Phi < DirectionMinPhi || Phi > DirectionMaxPhi)
                        {
                            continue;
                        }
						
						double CosTheta = uni();//only create x-rays with +ve z    0-->pi: uni()*2.0 - 1.0;
                        if(CosTheta < DirectionMinCosTheta || CosTheta > DirectionMaxCosTheta)
                        {
                            continue;
                        }
						
					
						double FluoEnergy = FluoData.PickFluorescenceEnergy();
						if(FluoEnergy < 0.0) //energy will be negative if the xray emission does not pass the fluorescence yield test
						{
							continue;
						}
						
						
						float AbsorptionLength = (-1.0f / AbsorbCoeff) * log(uni()); //how far xray travelled into material
                        
                        while(AbsorptionLength > PathLength)
                        {
                            //need to guarantee that the photon is absorbed within the crystal.
                            //This is a stupid way to do this... If the absorption probablility is very low there's a reasonable chance of hanging here for a very long time.
                            AbsorptionLength = (-1.0f / AbsorbCoeff) * log(uni());
                        }
                        //isn't this caculation equivalent to checking if the photon is absorbed in the first place?
                        //should check if it is quicker just to do this test. (it's not)
						
						
						float Depth = fabs( AbsorptionLength * Direction.z ); //how deep xray gets into crystal
						
						
						float AbsorbCoeffFluo = MuData.GetAbsorbCoeffDataPoint(EnergyToWavelength(FluoEnergy));
						
						
                        float ProbTransmit = exp( -1.0f * AbsorbCoeffFluo * fabs(Depth/CosTheta)  ); //Depth/CosTheta to determine path length back out of crystal (reflection geometry only)
                        
                        if(uni() > ProbTransmit) //if absorbed, continue to next photon
                        {
                            //this checks if the photon would actually be transmitted (pass back out through the material)
                            continue;
                        }
                        
                        Vector EmitDirection( acos(CosTheta), Phi, true, true);
                        
                        
                        Vector NewSource(Source.x - (Source.z/Direction.z)*Direction.x,
                                         Source.y - (Source.z/Direction.z)*Direction.y,
                                         0.0f);
                        
                        double CCDIntersectX = NewSource.x + ((CCDZ-NewSource.z)/(EmitDirection.z))*EmitDirection.x;
                        double CCDIntersectY = NewSource.y + ((CCDZ-NewSource.z)/(EmitDirection.z))*EmitDirection.y;
                        

                        FluoResults << CCDIntersectX << "\t" << CCDIntersectY << endl;
                        AdvFluoResults << NewSource.x << "\t" << NewSource.y << "\t" <<
                        EmitDirection.x << "\t" << EmitDirection.y << "\t" << EmitDirection.z << "\t" << FluoEnergy << endl;
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
    cout << "Proc " << ProcessorId << ": " <<  "Done! (" << int(time(NULL) - StartTime) << "s)" << endl;

    MPI_Barrier(MPI_COMM_WORLD);
    
    //cat DiffractResults_P0.txt DiffractResults_P1.txt DiffractResults_P2.txt DiffractResults_P3.txt > DiffractResults.txt
    if( ProcessorId == 0)
    {
        system(CreateConcatCommand("DiffractResults", ".txt").c_str());
        system(CreateConcatCommand("AdvDiffractResults", ".txt").c_str());
        system(CreateConcatCommand("FluoResults", ".txt").c_str());
        system(CreateConcatCommand("AdvFluoResults", ".txt").c_str());
		
		system("rm AdvDiffractResults_P*");
		system("rm DiffractResults_P*");
		system("rm AdvFluoResults_P*");
		system("rm FluoResults_P*");
		
    }

    MPI_Finalize();
    
    if(ProcessorId == 0)
    {
        cout << endl;
        cout << "/////////////////////////////////////////////" << endl;
        cout << "//Finished running Main Diffraction Program//" << endl;
        cout << "/////////////////////////////////////////////" << endl;
        cout << endl;        
    }
    
    return 0;
}

//**********************************************************
//Main program of Caulobacter crescentus cell cycle model
//
//For model details please refer to 
//
//A Stochastic Spatiotemporal Model of a Response-Regulator Network 
//in the Caulobacter crescentus Cell Cycle
//Fei Li , Kartik Subramanian, Minghan Chen, John J. Tyson, Yang Cao
//
//Physical Biology, special issue of 2015 q-bio 
//
//coded by: Fei Li
//felix@cs.vt.edu
//Dec 4, 2015
//***********************************************************

#include "mpi.h"

#include <cstdio>
using std::sprintf;

#include <sstream>
#include <string>
using std::string; 

#include <cstring>
using std::strcpy; 
using std::strcat;

#include <cstdlib>
using std::system;
using std::srand;

#include <ctime>
using std::time;

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include "PleC.h"

int main(int argc, char* argv[]){
 


    //*************************************
    //make directory for output
    //mkdir outDir
    //*************************************
    int stat=0;
    string mkoutDir = "mkdir -p ";
    string outDir = "./data/";
    mkoutDir.append(outDir);

    stat=system(mkoutDir.c_str());

    if (stat!=0) {
        cerr<<"Unable to creat directory: "<<mkoutDir<<endl;
        exit(1);
    }


    int   numtasks, rank, len;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //*************************************
    //Usage: kartik <filename>
    //*************************************
     std::ostringstream convert;   // stream used for the conversion
     convert << rank;         // insert the textual representation of 'Number' in the characters in the stream

     string filenum = convert.str();


    //*******************************************
    //the file to keep the population trajectory 
    //*******************************************
    string popDist, fireDist;
    popDist.append(outDir);
    popDist.append("TraFei");
    popDist.append(filenum);
    popDist.append(".txt");

    fireDist.append(outDir);
    fireDist.append("FireFei");
    fireDist.append(filenum);
    fireDist.append(".txt");

    PleC pc;

    //****************************************
    //seed the random variable 
    // Random Numbers In Scientific Computing: An Introduction by Katzgrabber 
    //****************************************
    srand(abs((time(NULL)*181)*((rank-83)*359)%104729));

    pc.InitStage();
    pc.initSpatial();
    pc.DNAorigination();
    pc.clearFire();

    time_t t1,t2;
    time(&t1);

    //********************************************
    //Fixed Length Stage for cell initiation
    //********************************************
    pc.cellFixed();
    pc.SSA(150.0, popDist, 1.0);
    pc.printFirings(fireDist);
    //*******************************************
    //Simulatio of wild type cells
    //*******************************************

    pc.clearFire();
    pc.cellGrowth(); 
    pc.SSA(30.0, popDist, 1.0);
    pc.printFirings(fireDist);

    pc.introDivJ();
    pc.SSA(20.0, popDist,1.0);
    pc.printFirings(fireDist);

    //********************************************
    //gene replication in 50 min
    //********************************************
    pc.DNAreplication();
    pc.clearPleC();
    pc.SSA(40.0, popDist,1.0);
    pc.printFirings(fireDist);

    pc.checkPoint3();
    pc.SSA(30.0, popDist,1.0);
    pc.printFirings(fireDist);

//    pc.compartmentization();
//    pc.SSA(30.0, popDist,1.0);

//pc.printFirings(fireDist);

    time(&t2);
    cout<<"Time Cost (hour)#"<<difftime(t2, t1)/3600.0<<endl;
    
    MPI_Finalize();

    return 0;
}

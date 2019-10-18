/*!
 * \file definitions.h
//author Yiming Li
//date 20170301
 * \brief The constants used throughout the project. 
 */
//#ifndef DEFINITIONS_H
//#define DEFINITIONS_H

#define BSZ 64			//!< max number of threads per block 
#define MEM 3872        //3872		11616//!< shared 
#define IMSZ 11			//!< not used
#define IMSZBIG 51		//!< maximum fitting window size
#define NK 128			//!< number of blocks to run in each kernel
#define pi 3.141592f	//!< ensure a consistent value for pi
#define NV_P 4			//!< number of fitting parameters for MLEfit (x,y,bg,I)
#define NV_PS 5			//!< number of fitting parameters for MLEFit_sigma (x,y,bg,I,Sigma)
#define NV_PZ 5			//!< not used (x,y,bg,I,z)
#define NV_PS2 6		//!< number of fitting parameters for MLEFit_sigmaxy (x,y,bg,I,Sx,Sy)
#define NV_PSP 5
//#define _DEBUG 1

#define TOLERANCE 1e-6f
#define INIT_ERR 1e13
#define INIT_LAMBDA 0.1f
#define SCALE_UP 10
#define SCALE_DOWN 0.1f
#define ACCEPTANCE 1.5f

#define BLOCK_MAX_SIZE 512


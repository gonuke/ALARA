/* $Id: Flux.C,v 1.14 2005-02-08 13:39:49 wilsonp Exp $ */
/* (Potential) File sections:
 * Service: constructors, destructors
 * Input: functions directly related to input of data 
 * xCheck: functions directly related to cross-checking the input
 *         against itself for consistency and completeness
 * Preproc: functions directly related to preprocessing of input
 *          prior to solution 
 * Solution: functions directly related to the solution of a (sub)problem
 * Utility: advanced member access such as searching and counting 
 */

#include <vector>

#include "Flux.h"
#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define SWAP_4(x) ( ((x) << 24) | (((x) << 8) & 0x00ff0000) | \
         (((x) >> 8) & 0x0000ff00) | ((x) >> 24) )
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_LONG(x) (*(unsigned *)&(x) = SWAP_4(*(unsigned *)&(x)))
#define FIX_FLOAT(x) FIX_LONG(x)

#include "Volume.h"
#include "Calc/VolFlux.h"

/***************************
 ********* Service *********
 **************************/

/** This constructor creates a blank list head when no arguments
    are given.  Otherwise, it sets the format, flux identifier, flux
    file name, scaling factor, and skip value, with arguments given in
    that order. */
Flux::Flux(int inFormat, int inSkip, double inScale, const std::string flxName, const std::string fName) :
  format(inFormat), skip(inSkip), scale(inScale), fluxName(flxName), fileName(fName)
{}

/** This assignmnet operator behaves similarly to the copy                        constructor. The correct implementation of this operator must ensure
    that previously allocated space is returned to the free store before          allocating new space into which to copy the object.  It does NOT              copy 'next'. */
Flux& Flux::operator=(const Flux& f) 
{
  if (this == &f)
    return *this;
  
  scale = f.scale;
  skip = f.skip;
  format = f.format;
  fluxName = f.fluxName;
  fileName = f.fileName;

  return *this;

}
/****************************
 *********** Input **********
 ***************************/

/***** get flux descriptions ******/
/* called by Input::read(...) */
/** It returns a pointer to the new object of class Flux which has just
    been created.  It does NOT read the actual flux information from              the file. */
Flux::Flux(istream& input)
{
  std::string type;
  
  input >> fluxName >> fileName >> scale >> skip >> type;

  switch(tolower(type[0]))
    {
    
    case 'r':
      // binary format (read from binary file)
      format = FLUX_R;
      break;
    
    case 'd':
      // default format (read from text file)
      format = FLUX_D;
      break;
    default:
      error(140,"Invalid flux type: %s", type.c_str());
    }

}

/****************************
 ********* Preproc **********
 ***************************/

/* cross-reference the fluxes with the intervals
 * read fluxes into interval member objects */
/* called by Input::preproc(...) */
/** The function expects a pointer to an object of class Volume which
    should be the head of the global interval list. */
void Flux::xRef(std::vector<Flux*> fluxList, Volume *volList)
{

  int numVols = volList->count();
  int numGrps = VolFlux::getNumGroups();
  double temp;

  VolFlux::setNumFluxes(fluxList.size());

  verbose(2,"Assigning %d fluxes to each interval",fluxList.size());

  // Dynamically Create Matrix
  double **FluxMatrix = new double*[numVols];
  double *MatrixStorage = new double[numVols*numGrps];

  for(int i = 0; i < numVols; i++)
    FluxMatrix[i] = &MatrixStorage[i*numGrps];

  /* for each flux definition */
  for ( std::vector<Flux*>::iterator ptr = fluxList.begin();
        ptr != fluxList.end();
        ptr++)
    {
      verbose(3,"Assigning flux %s",(*ptr)->fluxName.c_str());

      switch ((*ptr)->format)
        {
        case FLUX_D:
          {
            /* Default: Reads data from fluxin file */

            // Open Input File
            ifstream FluxData((*ptr)->fileName.c_str());

            // Skip appropriate number of Volumes
            if((*ptr)->skip > 0)
              for(int i = 0; i < (*ptr)->skip; i++)
                for(int j = 0; j < numGrps; j++)
                  FluxData >> temp;

            if(FluxData.eof())
              error(622,"Flux file %s does not contain enough data.",
                    (*ptr)->fileName.c_str());

            // Load data from FluxData
            for(int x = 0; x < numVols; x++)
              {
                for(int y = 0; y < numGrps; y++)
                  {
                    if(FluxData.eof())
                      error(622,"Flux file %s does not contain enough data.",
                            (*ptr)->fileName.c_str());

                    FluxData >> FluxMatrix[x][y];
                  }
              }
            FluxData.close();
            break;
          }

        case FLUX_R:
          {
            (*ptr)->readRTFLUX(MatrixStorage,numVols,numGrps);
            
            break;
          }
        };
      
      volList->storeMatrix(FluxMatrix,(*ptr)->scale);
    }
  
  delete[] FluxMatrix; // FIXME: internal arrays not deleted?

  verbose(3,"Assigned %d fluxes to each interval",fluxList.size());

}


void Flux::readRTFLUX(double *MatrixStorage,int numVols, int numGrps)
{
 
  FILE* binFile = fopen(fileName.c_str(),"rb");

  int f77_reclen; /// needed to accommodate strange F77 binary strucutre
  int readInt;
  float readFlt;
  char buffer[256];

  int grpLo, grpUp, grpHi;

  /// read file header
  fread((char*)&f77_reclen,SINT,1,binFile);
  debug(2,"readRTFLUX: f77_reclen: %d",f77_reclen);

  fread(buffer,1,24,binFile);
  fread((char*)&readInt,SINT,1,binFile);
  fread((char*)&f77_reclen,SINT,1,binFile);

  debug(2,"readRTFLUX: f77_reclen: %d",f77_reclen);

  /// read dimensions in file
  int ndim, ngrp, ninti, nintj, nintk, nblok;
  fread((char*)&f77_reclen,SINT,1,binFile);
  fread((char*)&ndim, SINT,1,binFile);
  fread((char*)&ngrp, SINT,1,binFile);
  fread((char*)&ninti,SINT,1,binFile);
  fread((char*)&nintj,SINT,1,binFile);
  fread((char*)&nintk,SINT,1,binFile);
  fread((char*)&readInt, SINT,1,binFile);
  fread((char*)&readFlt, SFLOAT,1,binFile);
  fread((char*)&readFlt,SFLOAT,1,binFile);
  fread((char*)&nblok, SINT,1,binFile);
  fread((char*)&f77_reclen,SINT,1,binFile);

  debug(2,"readRTFLUX: (ndim,ngrp,ninti,nintj,nintk,nblok) = (%d,%d,%d,%d,%d,%d)",
	ndim,ngrp,ninti,nintj,nintk,nblok);

  /// error checking
  if (ndim > 1)
    error(624,"RFLUX file: %s is 2- or 3-dimensional.  This feature currently only supports 1-D.",fileName.c_str());

  if (ngrp<numGrps)
    error(623,"RTFLUX file: %s does not contain enough data - not enough groups", fileName.c_str());

  if (ninti<(skip+numVols))
    error(623,"RTFLUX file: %s does not contain enough data - not enough intervals", fileName.c_str());

  /// read blocks (1-D)
  double* fluxIn = new double[ninti*ngrp];
  for (int blkNum=0;blkNum<nblok;blkNum++)
    {
      grpLo =   blkNum   * ((ngrp-1)/nblok + 1);
      grpUp = (blkNum+1) * ((ngrp-1)/nblok + 1)-1;
      grpHi = std::min(ngrp,grpUp);
      fread((char*)&f77_reclen,SINT,1,binFile);
      fread((char*)(fluxIn+grpLo*ninti),SDOUBLE,(grpHi-grpLo+1)*ninti,binFile);
      fread((char*)&f77_reclen,SINT,1,binFile);
    }
  
  debug(2,"readRTFLUX: reading %d groups in %d volumes, skipping %d entries", numGrps,numVols,skip);
  /// transpose data
  for (int gNum=0;gNum<numGrps;gNum++)
    for (int volNum=0;volNum<numVols;volNum++)
      {
	debug(3,"readRTFLUX: reading group #%d in volume #%d: %g", 
	      gNum, volNum, fluxIn[gNum*ninti+(volNum+skip)]);
	MatrixStorage[volNum*numGrps+gNum] = fluxIn[gNum*ninti+(volNum+skip)];
      }

  delete fluxIn;

  return;

}

/****************************
 ********* Utility **********
 ***************************/

/** Returns the 0-based ordinal number of the flux in the list.
    Special values (< 0) are returned for special cases, such as bad
    file names or unfound flux descriptions. */
int Flux::find(std::vector<Flux*> fluxList, const std::string srchFlux)
{
  for ( std::vector<Flux*>::iterator ptr = fluxList.begin();
        ptr != fluxList.end();
        ptr++)
    {
      if (!(*ptr)->fluxName.compare(srchFlux))
        {
          if ((*ptr)->checkFname())
            return ptr-fluxList.begin();
          else
            return FLUX_BAD_FNAME;
        }
    }
  
  return FLUX_NOT_FOUND;
}

int Flux::checkFname()
{
  FILE* textFile = fopen(fileName.c_str(),"r");

  if (textFile != NULL)
    {
      verbose(5,"Openned flux file %s.",fileName.c_str());
      fclose(textFile);
      return TRUE;
    }
  else
    {
      warning(340,"Unable to open flux file %s for flux %s.",fileName.c_str(),fluxName.c_str());
      return FALSE;
    }
 
}

/* $Id: Flux.h,v 1.13 2004-08-04 19:44:52 wilsonp Exp $ */
#include "alara.h"

#ifndef FLUX_H
#define FLUX_H

/* flux formats */
#define FLUX_HEAD 0
#define FLUX_D 1
#define FLUX_R 2

/* results of searching for flux */
#define FLUX_NOT_FOUND -1
#define FLUX_BAD_FNAME -2

/** \brief This class is invoked as a linked list and describes the 
 *         flux information of the problem.
 *   
 *  The single object of class Input has a list of Fluxes.  The first
 *  element in each list has type FLUX_HEAD (defined through the '
 *  format' member), and contains no problem data.
 */

class Flux
{
protected:
  /// This indicates the format of this flux and is based on the
  /// definitions given below.
  int format;
    
  /// The number of full interval entries to skip before the first one
  /// to be read.
  /** This makes it easy to either skip entries not modelled in the
      geometry, use the same file for a number of flux
      descriptions - subsequent descriptions need to skip the initial
      ones.  Again, note that this is a count of how many *full interval
      entries* - one entry has nGroup data values. */
  int skip;

  /// A scalar normalization for this flux description.
  /** This allows re-normalization of flux spectra either for physical
      or modelling reasons. */
  double scale;

  /// The descriptive name of this flux description.
  /** This will be used primarily to cross-reference this flux 
      description with the schedules. */
  std::string fluxName;
    
  /// The filename where these flux spectra should be read.
  std::string fileName;


public:
  /// Default constructor
  Flux(int inFormat=FLUX_HEAD, int inSkip=0, double Scale=0, 
       const std::string flxName="", const std::string fName="");

  /// Copy constructor
  Flux(const Flux&);

  /// Constructor to make a new Flux object from the input file.
  Flux(istream&);

  /// Inline destructor 
  /** Deletes storage for the flux identifier and the flux file name.
      The whole list is destroyed by deleting the next pointer. */

  /// Overloaded assignment operator
  Flux& operator=(const Flux&);

  /// Function to read RTFLUX binary FORTRAN output from DANTSYS, etc.
  void readRTFLUX(double*, int, int);

  /// Setup the problem to read the actual flux data into containers in each
  /// interval.
  static void xRef(std::vector<Flux*>, Volume*);
  
  /// Find a specific flux description based on the flux identifier give
  /// as the argument.
  static int find(std::vector<Flux*>, std::string);

  /// Check that the filename specified in the description is useful by
  /// openning and closing again.
  int checkFname();

};



#endif

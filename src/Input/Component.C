/* $Id: Component.C,v 1.18 2003-10-28 22:11:37 wilsonp Exp $ */
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


#include "Component.h"
#include "Mixture.h"

#include "Chains/Root.h"

ifstream Component::matLib;
std::map<std::string, Component::EleLibEntry> Component::eleLib;

/***************************
 ********* Service *********
 **************************/
/** This constructor creates a blank list head, when no arguments
    are given.  Otherwise, it sets the type, the name and the density.
    The 'next' element is initialized to NULL. */
Component::Component(int compType, std::string name, double dens, double volFrac) :
  type(compType),density(dens), volFraction(volFrac),compName(name)
{}


/** The assignment operator is similar to the copy constructor, but it
    uses an already allocated object on the left hand side.  The
    correct implementation of this operator must ensure that
    previously allocated space is returned to the free store before
    allocating new space into which to copy the object. Note that
    'next' is NOT copied, the left hand side object will continue to
    be part of the same list unless explicitly changed. */
Component& Component::operator=(const Component& comp)
{ 
  if (this != &comp)
    {
      type = comp.type;
      compName = comp.compName;
      density = comp.density;
      volFraction = comp.volFraction;
    }

  return *this;

}


/****************************
 *********** Input **********
 ***************************/

/** It expects an integer type, detemined by the calling function,
    and a reference to the input file's stream.  It returns
    pointer to the new object of class Component which has just been
    read. */
Component* Component::getComponent(int setType,istream &input, Mixture *mixPtr)
{

  std::string name;
  double dens=0, volFrac=0;

  input >> name;
  if (setType < COMP_SIM)
    input >> dens >> volFrac;
  if (setType == COMP_SIM)
    input >> volFrac;

  next = new Component(setType,name,dens,volFrac);
  memCheck(next,"Component::getComponent(...) : next");

  verbose(3,"type code: %d name: %s, density %g, volume fraction: %g",
          setType,name.c_str(),dens,volFrac);

  mixPtr->incrVolFrac(volFrac);

  return next;

}

void Component::getMatLib(istream& input)
{
  std::string fname;
  input >> fname;
  matLib.open(searchNonXSPath(fname.c_str()));

  if (matLib == 0)
    error(110,"Unable to open material library: %s",fname.c_str());

  verbose(2,"Openned material library %s",searchNonXSPath(fname.c_str()));
}

void Component::getEleLib(istream& input)
{
  std::string fname;
  input >> fname;
  ifstream eleLibFile;

  eleLibFile.open(searchNonXSPath(fname.c_str()));

  if (eleLibFile == 0)
    error(110,"Unable to open element library: %s",fname.c_str());

  verbose(2,"Openned element library %s",searchNonXSPath(fname.c_str()));

  std::string eleKey;
  Component::EleLibEntry blank_element;
  std::pair<std::string,double> isoInfo;

  while (!eleLibFile.eof())
    {
      clearComment(eleLibFile);
      eleLibFile >> eleKey;
      eleLib[eleKey] = blank_element;
      eleLibFile >>  eleLib[eleKey].A >> eleLib[eleKey].Z 
                 >> eleLib[eleKey].density >> eleLib[eleKey].numIsos;
      int numIsos = eleLib[eleKey].numIsos;
      while (numIsos-->0)
        {
          clearComment(eleLibFile);
          eleLibFile >> isoInfo.first >> isoInfo.second;
          eleLib[eleKey].isoList.push_back(isoInfo);
        }
    }
 
  eleLibFile.close();
 
}


/****************************
 ********* Preproc **********
 ***************************/

/** It expects a list of components, through a pointer to the head
    of the list.  Rather than deleting the current component, it is
    changed to a copy of the first object, and the others in the list
    are inserted as new objects.  This is used to replace all 'similar'
    components before a mixture is expanded. */
Component* Component::replaceSim(Component *newCompList)
{
  Component *ptr = this;
  Component *saveNext = ptr->next;
  Component *newComp = newCompList->next;
  double scale = volFraction;

  /* change this component */
  *ptr = *newComp;
  ptr->volFraction *= scale;

  /* insert all subsequent components from this mixture */
  while (newComp->next != NULL)
    {
      newComp = newComp->next;
      ptr->next = new Component(*newComp);
      memCheck(ptr->next,"Component::replaceSim(...) : ptr->next");
      ptr = ptr->next;
      ptr->volFraction *= scale;
    }

  ptr->next = saveNext;

  return ptr;
  
}

/** It is always called through the head of the Component list for a
    given mixture.  For cross-referencing, a pointer to that mixture
    is expected as an argument. */
Root* Component::expand(Mixture *mix)
{
  Component* ptr = this;
  Root* rootList = new Root;
  Root* compRootList;
  memCheck(rootList,"Component::expand(...) : rootList");

  while (ptr->next != NULL)
    {
      ptr = ptr->next;
      
      switch (ptr->type)
	{
	case COMP_MAT:
	  compRootList = ptr->expandMat(mix);
	  rootList = rootList->merge(compRootList);
	  verbose(6,"Merged material %s into rootList for mixture",
              ptr->compName.c_str());
	  delete compRootList;
	  break;
	case COMP_ELE:
	case TARGET_ELE:
	  compRootList = ptr->expandEle(mix,ptr);
	  rootList = rootList->merge(compRootList);
	  verbose(6,"Merged element %s into rootList for mixture",
              ptr->compName.c_str());
	  delete compRootList;
	  break;
//	case COMP_ISO:
	case TARGET_ISO:
	  Root* newRoot = new Root(ptr->compName.c_str(),ptr->density,mix,ptr);
	  memCheck(newRoot,"Component::expand(...) : newRoot");
	  rootList = rootList->merge(newRoot);
	  verbose(6,"Merged isotope %s into rootList for mixture",
              ptr->compName.c_str());
	  delete newRoot;
	  break;
	}
    }
  return rootList;

}

/** For cross-referencing, it expects a pointer to the mixture and
    component which contain this particular element - it does not
    automatically use the 'this' component pointer, since it can be
    called through a temporary object (such as might be created in
    expandMat() while expanding a material).  It returns a pointer to
    an object of class Root, which serves as the head of a list. */
Root* Component::expandEle(Mixture* mix, Component* comp)
{
  Root *rootList = new Root;
  memCheck(rootList,"Component::expandEle(...) : rootList");

  verbose(4,"Expanding element %s",compName.c_str());

  if (eleLib.count(compName) != 0)
    {
      Component::EleLibEntry ele = eleLib[compName];
      
      if (density >= 0)
        density *= ele.density;
      else
        density = -density;
      
      double Ndensity = volFraction * density * AVAGADRO/ele.A;
      mix->incrTotalDensity(density*volFraction);
      
      /* if element is found, add a new root for each isotope */
      verbose(5,"Found element %s with %d isotopes in element library",
              compName.c_str(), ele.numIsos);
      
      for (std::vector<std::pair<std::string,double> >::iterator iso = ele.isoList.begin();
           iso != ele.isoList.end();
           iso++)
        {
          std::string isoName = compName.substr(compName.find(':')) + '-' + iso->first;
          Root* newRoot = new Root(isoName.c_str(),iso->second*Ndensity/100.0,mix,comp);
          memCheck(newRoot,"Component::expandEle(...) : newRoot");
          rootList = rootList->merge(newRoot);
          verbose(6,"Merged isotope %s into rootList for element %s",
                  isoName.c_str(),compName.c_str());
          debug(5,"Accounted for isotope %s in Root List",isoName.c_str());
          delete newRoot;
        }
    }
  else
    error(310,"Could not find element %s in element library.",compName.c_str());
  
  return rootList;

}

/** For cross-referencing, it expects a pointer to the mixture and
    component which contain this particular element - it uses the
    'this' component pointer since it will never be ambiguous.  It
    returns a pointer to an object of class Root, which serves as the
    head of a list. */
Root* Component::expandMat(Mixture* mix)
{
  Root *rootList = new Root;
  memCheck(rootList,"Component::expandMat(...) : rootList");
  Component *element;
  int numEles, eleZ;
  char testName[64],eleName[64];
  double eleDens, matDens;

  /* rewind the material library */
  matLib.seekg(0L,ios::beg);

  verbose(4,"Expanding material %s",compName.c_str());

  /* search for this material */
  clearComment(matLib);
  matLib >> testName >> matDens >> numEles;
  while (strcmp(testName,compName.c_str()) && !matLib.eof())
    {
      verbose(5,"Skipping material %s in material library.",testName);
      while (numEles-->0)
	{
	  clearComment(matLib);
	  matLib >> eleName >> eleDens >> eleZ;
	}
      clearComment(matLib);
      matLib >> testName >> matDens >> numEles;
    }

  if (!matLib.eof())
    {
      /* scale relative density by material density from lib */
      density *= matDens;

      verbose(5,"Found material %s in material library.",testName);
      /* if material found, read list of elements,
       * supplementing the root list for each one */
      while (numEles-->0)
	{
	  clearComment(matLib);
	  matLib >> eleName >> eleDens  >> eleZ;
	  eleDens *= -density*volFraction/100.0;
	  element = new Component(COMP_ELE,eleName,eleDens);
	  memCheck(element,"Component::expandMat(...) : element");
	  Root *elementRootList = element->expandEle(mix,this);
	  rootList = rootList->merge(elementRootList);
	  delete elementRootList;
	  verbose(6,"Merged element %s into rootList for material %s",
              eleName,compName.c_str());
	  delete element;
	}
    }
  else
    error(311,"Could not find material %s in material library.",
          compName.c_str());
  
  return rootList;
}

/****************************
 ********* Utility **********
 ***************************/

int Component::getCompNum(Component *compPtr)
{
  int compNum = -1;
  Component *ptr = this;

  while (ptr != NULL && ptr != compPtr)
    {
      ptr = ptr->next;
      compNum++;
    }

  return compNum;
}

/** It is used to primarily to search for Components of type 'similar'.
    It expects an integer argument giving the component type of
    interest. */
Component* Component::exists(int srchType)
{
 
  Component* ptr = this;

  while (ptr->next != NULL)
    {
      ptr = ptr->next;
      if (ptr->type == srchType)
	return ptr;
    }

  /* return the NULL to denote the end of this component list */
  return ptr->next;
}


  

//
//  EssentialPatchBC.hpp
//  lifev-heart
//
//  Created by Thomas Kummer on 27.09.17.
//  Copyright © 2017 Thomas Kummer. All rights reserved.
//

#ifndef EssentialPatchBC_hpp
#define EssentialPatchBC_hpp

#include <stdio.h>
#include <lifev/em/examples/example_EMHeart/GenericFactory.hpp>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_MpiComm.h"

namespace LifeV
{

class EssentialPatchBC
{
public:
    
    typedef VectorEpetra                                    vector_Type;
    typedef boost::shared_ptr<vector_Type>                  vectorPtr_Type;
    
    typedef BCVector                                        bcVector_Type;
    typedef boost::shared_ptr<bcVector_Type>                bcVectorPtr_Type;
    
    typedef EssentialPatchBC                                super;

    
    EssentialPatchBC(){}
    ~EssentialPatchBC(){}
    //
    virtual void setup(const GetPot& dataFile, const std::string& name) //so here we read from a file and assign some values
    {
    	
        // Patch name
        m_Name = name;
		
	//std::cout << "This is variable m_Name: " << m_Name << std::endl;
		        
        // Epicardium flag
        m_PrevFlag = dataFile ( ("solid/boundary_conditions/" + m_Name + "/flag").c_str(), 0.0 ); //what does zero stand for?
        
        // Patch motion direction   //what do we mean by patch direction? with respect to which coordinate system?
        for ( UInt j (0); j < 3; ++j ) //j (0) is just another way to set j = 0
        {
            m_patchDirection[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/direction").c_str(), 0.0, j );
        }
        m_patchDirection.normalize();
        
        // Boundary condition components //I don't get this one here
        UInt componentSize = dataFile.vector_variable_size ( ("solid/boundary_conditions/" + m_Name + "/component").c_str() );
        for ( UInt j (0); j < componentSize; ++j )
        {
            m_patchComponent.push_back( dataFile ( ("solid/boundary_conditions/" + m_Name + "/component").c_str(), 0.0, j ) );
        }
        

        // Patch peak displacement
        m_patchDisplacement = dataFile ( ("solid/boundary_conditions/" + m_Name + "/displacement").c_str(), 1.0 );
        
        // Temporal activation parameter
        m_tmax = dataFile ( "solid/patches/tmax", 0. );
        m_tduration = dataFile ( "solid/patches/tduration", 0. );
    }
    
    void createPatchArea (EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const int& newFlag, const Real& time)
    {
    	auto p2FeSpace = solver.electroSolverPtr()->feSpacePtr();
    	auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
    	FESpace<RegionMesh<LinearTetra>, MapEpetra > p1FESpace (p2FeSpace->mesh(), "P1", 1, p2FeSpace->mesh()->comm());

    	VectorEpetra p1ScalarFieldFaces (p1FESpace.map());
    	p1ScalarFieldFaces *= 0.0;
    	Int p1ScalarFieldFacesDof = p1ScalarFieldFaces.epetraVector().MyLength();
    	//int globalIdArray[p1ScalarFieldFacesDof];

    	 m_patchFlag = newFlag;
    	 const auto& mesh = solver.localMeshPtr();

    	 int globalIdArray[p1ScalarFieldFacesDof];
    	 p1ScalarFieldFaces.blockMap().MyGlobalElements(globalIdArray);


    	 for (int j(0); j < mesh->numBoundaryFacets(); j++) //returns number of boundary facets
    	 {
    	       auto& face = mesh->boundaryFacet(j);
    	       auto faceFlag = face.markerID();
		
		//std::cout << "This is faceFlag in createPatchArea before changing: " << faceFlag << std::endl;
		//std::cout << "This is value of m_PrevFlag (should be 464): " << m_PrevFlag << std::endl;


    	       if (faceFlag == m_PrevFlag)
    	       {
    	            int numPointsOnFace(0);

    	            for (int k(0); k < 3; ++k)
    	            {

    	               	ID pointGlobalId = face.point(k).id();
    	               	auto coord = face.point(k).coordinates();
    	              // 	auto pointInPatch = nodeOnPatch

    	               	auto pointInPatch = nodeOnPatchCurrent(coord, time);

    	               	if(pointInPatch == true)
    	                {
    	               		++numPointsOnFace;
    	                  	for(int n = 0; n < p1ScalarFieldFacesDof; n++)
    	                  	{
    	                  		if(pointGlobalId == globalIdArray[n])
    	                  		{
    		                  		p1ScalarFieldFaces[pointGlobalId] = 1.0;
	
    	                  		}
    	                  	}
    	               	}

    	              }


    	 //if (numPointsOnFace > 2) // if there are more than two points on face we execute the if statement; not completly sure here

    	 	      if(numPointsOnFace >= 1)
    	          {
    	                     face.setMarkerID(m_patchFlag);
				auto faceFlagChanged = face.markerID();
				//std::cout << "This is changed faceFlag in createPatchArea: " << faceFlagChanged << std::endl;

    	          }
    	  }
    	 }

    	 m_patchFacesLocationPtr.reset(new vector_Type (p2FeSpace->map() ));
    	 *m_patchFacesLocationPtr = p2FeSpace->feToFEInterpolate(p1FESpace, p1ScalarFieldFaces);

    	 m_patchLocationPtr.reset (new vector_Type (p2FeSpace->map() ));
    }



    void applyBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const GetPot& dataFile)
    {
        auto dFeSpace = solver.structuralOperatorPtr() -> dispFESpacePtr(); 
      
        m_dispPtr = solver.structuralOperatorPtr()->displacementPtr(); 
		
        m_patchDispPtr = directionalVectorField(solver,dFeSpace, m_patchDirection, 1e-10, 0.0);
	
	*m_patchDispPtr *= 0.0; 
	
        m_patchDispBCPtr = bcVectorPtr_Type( new bcVector_Type( *m_patchDispPtr, dFeSpace -> dof().numTotalDof(), 1 ) );
        solver.bcInterfacePtr() -> handler()->addBC (m_Name, m_patchFlag,  Essential, Component, *m_patchDispBCPtr, m_patchComponent);
       
    }
    
    
    void modifyPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time, int& PatchFlag)
    {
	
	const int currentPatchFlag = PatchFlag;

        auto dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
        
        modifyPatchArea(solver, currentPatchFlag, time);

        Real currentPatchDisp = activationFunction(time) + 1e-3;
        if ( 0 == solver.comm()->MyPID() ) std::cout << "\nEssentialPatchBC: " << m_Name << " displaced by " << currentPatchDisp << " cm";

        m_patchDispPtr = directionalVectorField(solver,dFeSpace, m_patchDirection, currentPatchDisp, time);

        m_patchDispBCPtr.reset( new bcVector_Type( *m_patchDispPtr, dFeSpace->dof().numTotalDof(), 1 ) );
  	solver.bcInterfacePtr()->handler()->modifyBC(currentPatchFlag, *m_patchDispBCPtr); //This was the version how it worked
    }
    
    
    vector_Type patchDisplacement(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
        auto dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
        vector_Type localPatchDisplacement ( dFeSpace->map(), Repeated );
        localPatchDisplacement *= 0.0;

        auto nCompLocalDof = localPatchDisplacement.epetraVector().MyLength() / 3;

        for (int j (0); j < nCompLocalDof; ++j) 
        {
            UInt iGID = m_patchDispPtr->blockMap().GID (j);
            UInt jGID = m_patchDispPtr->blockMap().GID (j + nCompLocalDof);
            UInt kGID = m_patchDispPtr->blockMap().GID (j + 2 * nCompLocalDof);

            localPatchDisplacement[iGID] = (*m_dispPtr)[iGID]; // * (*m_patchLocationPtr)[iGID];
            localPatchDisplacement[jGID] = (*m_dispPtr)[jGID]; // * (*m_patchLocationPtr)[iGID];
            localPatchDisplacement[kGID] = (*m_dispPtr)[kGID]; // * (*m_patchLocationPtr)[iGID];
        }

        return localPatchDisplacement;
    }

    vector_Type displayDirectionalVectorField(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {

    	Vector3D current_point_on_plane;
    	Real distance;
    	Vector3D normalVector;
    	Vector3D startingPoint;
    	Vector3D direction = normalVector;
    	direction.normalize();


    	startingPoint[0] = -3.76487;
    	startingPoint[1] = -10.6687;
    	startingPoint[2] = -0.36572;

    	normalVector[0] = 0.665647;
    	normalVector[1] = 0.695607;
    	normalVector[2] = -0.270367;

    	//first we want to set up the initial vector
    	auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
    	  	const auto& meshFull = solver.fullMeshPtr();
    	    	FESpace<RegionMesh<LinearTetra>, MapEpetra > p1dFESpace (p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm());

    	    	    	VectorEpetra p1PositionVector (p1dFESpace.map());
    	    	    	p1PositionVector *= 0.0;

    	    	    	Int p1nLocalPoints = p1PositionVector.epetraVector().MyLength() / 3;

    	    	    	//std::cout << "This is length of positionVector: " << p1PositionVector.epetraVector().MyLength() << std::endl;

    	    	    	for (int j (0); j < p1nLocalPoints; j++)
    	    	    	{
    	    	    	     UInt iGID = p1PositionVector.blockMap().GID (j);
    	    	    	     UInt jGID = p1PositionVector.blockMap().GID (j + p1nLocalPoints);
    	    	    	     UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nLocalPoints);


    	    	    	     //Vector3D coord = meshFull->point(iGID).coordinates();

    	    	    	     p1PositionVector[iGID] = meshFull->point(iGID).x();
    	    	    	     p1PositionVector[jGID] = meshFull->point(iGID).y();
    	    	    	     p1PositionVector[kGID] = meshFull->point(iGID).z();



    	    	    	}

    	    	    	VectorEpetra p2PositionVector ( m_dispPtr->map() );
    	    	    	p2PositionVector = p2dFeSpace->feToFEInterpolate(p1dFESpace, p1PositionVector);

    	    	    	//now we have the positionvector


    	    	    	if(time == 0.0)
    	    	    	{
    	    	    		m_currentPositionVector = vectorPtr_Type (new VectorEpetra( p1dFESpace.map(), Repeated ));
    	    	    	    //*m_currentPositionVector = p2PositionVector;
    	    	    	}


    	    	    	vectorPtr_Type p2PatchDisplacement (new VectorEpetra( p1dFESpace.map(), Repeated ));
    	    	    	auto nCompLocalDof = p2PatchDisplacement->epetraVector().MyLength() / 3;




    	    	    	if(normalVector[2] != 0.0)
    	    	    		        {
    	    	    		        	//In thoughts we set cooridnates x and y equal to zero and solve for z coordinate and store it in current_point_on_plane[0]
    	    	    		        	//here we just added the max value of distance vector (3.5155), let's see how it works
    	    	    		        	current_point_on_plane[2] = (normalVector[0]*startingPoint[0] + normalVector[1]*startingPoint[1] + normalVector[2]*startingPoint[2] +activationFunction(time))/normalVector[2];
    	    	    		        	current_point_on_plane[1] = 0;
    	    	    		        	current_point_on_plane[0] = 0;

    	    	    		        	//std::cout << "This is coordinate of current point on plane" << current_point_on_plane[2] << std::endl;

    	    	    	 	        }


    	    	    	 for (int j (0); j < nCompLocalDof; ++j)
    	    	    		        {
    	    	    		                    // Get coordinates

    	    	    		                    UInt iGID = p2PatchDisplacement->blockMap().GID (j);
    	    	    		                    UInt jGID = p2PatchDisplacement->blockMap().GID (j + nCompLocalDof);
    	    	    		                    UInt kGID = p2PatchDisplacement->blockMap().GID (j + 2 * nCompLocalDof);


    	    	    		                    Vector3D coordinates;

    	    	    		                    coordinates(0) = p2PositionVector[iGID];
    	    	    		                    coordinates(1) = p2PositionVector[jGID];
    	    	    		                    coordinates(2) = p2PositionVector[kGID];

    	    	    		                    Vector3D QP; //define here the vector that goes from Q (point on plane) to point P

    	    	    		                    QP = coordinates - current_point_on_plane;

    	    	    		                    if(QP.dot(direction) <= 0) //here i have change to > 0
    	    	    		                    {
    	    	    		                    			distance = 1.0;
    	    	    		                               	//distance = abs(QP.dot(direction));
    	    	    		                    }
    	    	    		                    else
    	    	    		                      {
    	    	    		                              	distance = 0.0;
    	    	    		                      }

    	    	    		                    Vector3D displacement_vector;
    	    	    		                    displacement_vector[0] = distance*direction[0];
    	    	    		                    displacement_vector[1] = distance*direction[1];
    	    	    		                    displacement_vector[2] = distance*direction[2];


    	    	    		                    (*p2PatchDisplacement)[iGID] = displacement_vector[0];
    	    	    		                    (*p2PatchDisplacement)[jGID] = displacement_vector[1];
    	    	    		                    (*p2PatchDisplacement)[kGID] = displacement_vector[2];

    	    	    		                    	                    (*m_currentPositionVector)[iGID] = (*m_currentPositionVector)[iGID] + (*p2PatchDisplacement)[iGID];
    	    	    		                    	                    (*m_currentPositionVector)[jGID] = (*m_currentPositionVector)[jGID] + (*p2PatchDisplacement)[jGID];
    	    	    		                    	                    (*m_currentPositionVector)[kGID] = (*m_currentPositionVector)[kGID] + (*p2PatchDisplacement)[kGID];

    	    	    		        }

    	    	    	 return *p2PatchDisplacement;
    }

    vector_Type patchVectorField(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {


    	auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
    	const auto& meshFull = solver.fullMeshPtr();
    	FESpace<RegionMesh<LinearTetra>, MapEpetra > p1dFESpace (p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm());

    	VectorEpetra p1VectorField (p1dFESpace.map());
    	p1VectorField *= 0.0;

    	Int p1nLocalPoints = p1VectorField.epetraVector().MyLength() / 3;

    	for (int j (0); j < p1nLocalPoints; j++)
    	{
    	     UInt iGID = p1VectorField.blockMap().GID (j);
    	     UInt jGID = p1VectorField.blockMap().GID (j + p1nLocalPoints);
    	     UInt kGID = p1VectorField.blockMap().GID (j + 2 * p1nLocalPoints);


    	     Vector3D coord = meshFull->point(iGID).coordinates();
    	     bool pointInPatch = nodeOnPatch(coord, time);

    	     if(pointInPatch == true)
    	     {
    	    	 p1VectorField[iGID] = 1.0;
    	    	 p1VectorField[jGID] = 0.0; //1.0;
    	    	 p1VectorField[kGID] = 0.0; //1.0;

    	     }



    	}

    	return p1VectorField;


   }

    vector_Type patchLocation()
    {
        return *m_patchLocationPtr;
    }
    
    vector_Type patchFacesLocation()
    {
    	return *m_patchFacesLocationPtr;
    }

protected:
    

    virtual vectorPtr_Type directionalVectorField (EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time)
    {

        vectorPtr_Type vectorField (new VectorEpetra( dFeSpace->map(), Repeated )); 
        auto nCompLocalDof = vectorField->epetraVector().MyLength() / 3; 

        direction.normalize(); //here the input direction vector gets normalized

        direction *= disp; //then the direction vector gets multiplied by the displacement
        
        for (int j (0); j < nCompLocalDof; ++j) //What happens in this for loop?
        {
        	//does GID stand for Group identification? Or I think now that it stands for global index

        	//what I don't get yet is the one with GID

            UInt iGID = vectorField->blockMap().GID (j); //UInt is just unsigned integer 32 bit
            UInt jGID = vectorField->blockMap().GID (j + nCompLocalDof);
            UInt kGID = vectorField->blockMap().GID (j + 2 * nCompLocalDof);


            //so here we dereference our pointer but somehow I don't get it
            (*vectorField)[iGID] = direction[0]; // direction is just a 3D vector with direction[0] is x coordinate of vector; are we telling here the different nodes how much they need to displace? //here I should test that with my own vector programm and see how it goes
            (*vectorField)[jGID] = direction[1];
            (*vectorField)[kGID] = direction[2];
        }

        return vectorField;
    }


    virtual void modifyPatchArea(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const int& newFlag, const Real& time){}

    virtual vector_Type p2PositionVectorInitial(const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> p2dFeSpace, EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver) const
    {
    	const auto& meshFull = solver.fullMeshPtr();
    	FESpace<RegionMesh<LinearTetra>, MapEpetra > p1dFESpace (p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm());

    	VectorEpetra p1PositionVector (p1dFESpace.map());
    	p1PositionVector *= 0.0;

    	Int p1nLocalPoints = p1PositionVector.epetraVector().MyLength() / 3;


    	for (int j (0); j < p1nLocalPoints; j++)
    	{
    	     UInt iGID = p1PositionVector.blockMap().GID (j);
    	     UInt jGID = p1PositionVector.blockMap().GID (j + p1nLocalPoints);
    	     UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nLocalPoints);


    	     p1PositionVector[iGID] = meshFull->point(iGID).x();
    	     p1PositionVector[jGID] = meshFull->point(iGID).y();
    	     p1PositionVector[kGID] = meshFull->point(iGID).z();



   	}

    	    VectorEpetra p2PositionVector ( m_dispPtr->map() );
    	    p2PositionVector = p2dFeSpace->feToFEInterpolate(p1dFESpace, p1PositionVector);

    	   return p2PositionVector;
    }




    virtual Real activationFunction (const Real& time) const
    {
        Real timeInPeriod = fmod(time - m_tmax + 0.5*m_tduration, 800.);
        bool inPeriod ( timeInPeriod < m_tduration && timeInPeriod > 0);
        Real sinusSquared = std::pow( std::sin(timeInPeriod * PI / m_tduration) , 2 ) * m_patchDisplacement;
        return ( inPeriod ? sinusSquared : 0 );
    }
    

    virtual const bool nodeOnPatch(const Vector3D& coord,const Real& time){} ; //this is my function
    //virtual const bool nodeOnPatch(const Vector3D& coord) const = 0;

    virtual const bool nodeOnPatchCurrent(const Vector3D& coord,const Real& time) {};
	
    virtual void initialisePositionVector(const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver){};

    
    std::string m_Name;
    unsigned int m_PrevFlag;
    unsigned int m_patchFlag;
    unsigned int m_flagIncreaserOne = 900;
    unsigned int m_flagIncreaserTwo = 901;   

 
    Vector3D m_patchDirection;
    Vector3D normal_vector;
    Vector3D starting_point;
    Real m_patchDisplacement;
    
    std::vector<ID> m_patchComponent;

    vectorPtr_Type m_dispPtr;

    vectorPtr_Type m_patchLocationPtr;
    vectorPtr_Type m_patchFacesLocationPtr;
    vectorPtr_Type m_currentPositionVector;

    vectorPtr_Type m_patchDispPtr;
    bcVectorPtr_Type m_patchDispBCPtr;


    Real m_tmax;
    Real m_tduration;
    
};

    

class EssentialPatchBCHandler
{
public:


    EssentialPatchBCHandler(const std::string& patchListName, const GetPot& dataFile) :
        m_patchListName ("solid/boundary_conditions/" + patchListName),
        m_dataFile (dataFile),
       m_patchNumber (( m_dataFile.vector_variable_size(m_patchListName.c_str()) ))
	{}
    
    ~EssentialPatchBCHandler(){}

    void addPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {

        for ( UInt i (0) ; i < m_patchNumber ; ++i )
        {
            const std::string patchName = m_dataFile ( m_patchListName.c_str(), " ", i ); //type is for example EssentialPatchBCCircular
            const std::string patchType = m_dataFile ( ("solid/boundary_conditions/" + patchName + "/type").c_str(), "EssentialPatchBCMovingPlane" );

            m_patchBCPtrVec.push_back(CREATE(EssentialPatchBC, patchType));
            m_patchBCPtrVec[i]->setup(m_dataFile, patchName);
            m_patchBCPtrVec[i]->createPatchArea(solver, 900 + i, time);
        }

        if ( solver.comm()->MyPID() == 0 ) std::cout << "\nEssentialPatchBCHandler: " << __FUNCTION__ << " - done" << std::endl;
    }

    void applyPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {

        m_patchDisplacementVecSumPtr = vectorPtr_Type (new VectorEpetra( solver.structuralOperatorPtr()->dispFESpacePtr()->map(), Repeated ));
        m_patchVecSumPtr = vectorPtr_Type (new VectorEpetra( solver.structuralOperatorPtr()->dispFESpacePtr()->map(), Repeated ));
        m_patchLocationScalarSumPtr = vectorPtr_Type (new VectorEpetra( solver.electroSolverPtr()->potentialPtr()->map(), Repeated ));
        m_patchFacesLocationScalarSumPtr = vectorPtr_Type (new VectorEpetra( solver.electroSolverPtr()->potentialPtr()->map(), Repeated ));
        m_directionVecFieldPtr =  vectorPtr_Type (new VectorEpetra( solver.structuralOperatorPtr()->dispFESpacePtr()->map(), Repeated ));

        for (auto& patch : m_patchBCPtrVec)
        {
            patch->applyBC(solver, m_dataFile);
        }

        updatePatchDisplacementSum(solver);
        updatePatchLocationSum(solver);
        updatePatchFacesLocationSum(solver);
        updatePatchVec(solver, 0.0);
        updateDispDirectionalVec(solver, 0.0);
        
        if ( solver.comm()->MyPID() == 0 ) std::cout << "\nEssentialPatchBCHandler: " << __FUNCTION__ << " - done" << std::endl;
    }

    void modifyPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {
    	int PatchFlag = 899;
	//std::cout << PatchFlag << std::endl;
        for (auto& patch : m_patchBCPtrVec)
        {
        	PatchFlag += 1;
            patch->modifyPatchBC(solver, time, PatchFlag);
        }

        updatePatchDisplacementSum(solver);
        updatePatchLocationSum(solver);
        updatePatchFacesLocationSum(solver);
        updatePatchVec(solver, time);
        updateDispDirectionalVec(solver, time);
        
        if ( solver.comm()->MyPID() == 0 ) std::cout << "\nEssentialPatchBCHandler: " << __FUNCTION__ << " - done" << std::endl;
    }

///////// THE FOLLOWING THREE FUNCTIONS WERE JUST USE TO WRITE COORDINATES TO FILES; WAS USED TO DEBUG; NOT NEEDED ANYMORE NOW
/*
    void getcoordinates(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
           	auto p2dFeSpace = solver.structuralOperatorPtr() -> dispFESpacePtr();
           	FESpace<RegionMesh<LinearTetra> , MapEpetra > p1dFESpace ( p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm() );

           	VectorEpetra p1PositionVector (p1dFESpace.map());

           	        // Fill P1 vector with mesh values
           	        Int p1nCompLocalDof = p1PositionVector.epetraVector().MyLength() / 3;

           	        //This works to get number of processors
           	        //std::cout << "THIS IS NUMBER OF PROCESSORS: " << p2dFeSpace->mesh()->comm()->NumProc() << std::endl;
           	        	//	std::cout << "" << std::endl;

           	        //const char* path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinates.dat";

           	        for (int j (0); j < p1nCompLocalDof; j++)
           	        {
           	            UInt iGID = p1PositionVector.blockMap().GID (j);
           	            UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
           	            UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);

           	            ID local_id = p2dFeSpace->mesh()->point(j).localId();

           	            auto currentprocessor = p2dFeSpace->mesh()->comm()->MyPID();
           	            std::ostringstream oss;//this is to convert int to string
           	            oss << currentprocessor;

           	            std::string path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinatefiles/coordinates_" + oss.str() + ".dat";
           	            std::ofstream writer(path.c_str(), std::ios_base::app);
           	            if(writer.is_open()==false)
           	            {
           	            	std::cout << "error occured while opening the file" << std::endl;
           	            }
           	            writer << p2dFeSpace->mesh()->point(local_id).x() << "\t\t" << p2dFeSpace->mesh()->point(local_id).y() << "\t\t" << p2dFeSpace->mesh()->point(local_id).z() << std::endl;
           	            writer.close();

           	        }
    }


    void getp1PositionVector(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
    	//here in the first part we just create the p1postionvector; this is identical to part in p2position vector
    	auto p2dFeSpace = solver.structuralOperatorPtr() -> dispFESpacePtr();
    	FESpace<RegionMesh<LinearTetra> , MapEpetra > p1dFESpace ( p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm() );

    	// Create P1 VectorEpetra
    	VectorEpetra p1PositionVector (p1dFESpace.map());

    	// MyLength returns the local vector length on the calling processor
    	Int p1nCompLocalDof = p1PositionVector.epetraVector().MyLength() / 3;

    	for (int j (0); j < p1nCompLocalDof; j++)
    	{
    		UInt iGID = p1PositionVector.blockMap().GID (j);
    		UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
    		UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);

    		ID local_id = p2dFeSpace->mesh()->point(j).localId(); //Here we get the localId of the point


    		p1PositionVector[iGID] = p2dFeSpace->mesh()->point (local_id).x();
    		p1PositionVector[jGID] = p2dFeSpace->mesh()->point (local_id).y();
    		p1PositionVector[kGID] = p2dFeSpace->mesh()->point (local_id).z();

    	}

    	//Now we have created the p1positionvector and now we want to print it out --> similar to getcoordinates
    	for(int j(0); j < p1nCompLocalDof; j++)
    	{
    		UInt iGID = p1PositionVector.blockMap().GID (j);
    		UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
    		UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);

    		auto currentprocessor = p2dFeSpace->mesh()->comm()->MyPID();
    		std::ostringstream oss;//this is to convert int to string
    		oss << currentprocessor;

       		std::string path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/positionvectorfiles/coordinates_" + oss.str() + ".dat";
    		std::ofstream writer(path.c_str(), std::ios_base::app);
    		if(writer.is_open()==false)
       		{
    			std::cout << "error occured while opening the file" << std::endl;
    		}
       		writer << p1PositionVector[iGID] << "\t\t" << p1PositionVector[jGID] << "\t\t" << p1PositionVector[kGID] << std::endl;
      		writer.close();

    	}

    }


    void getcoordinates_createPatchArea(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const int& newFlag, const Real& time)
    {

    	unsigned int m_patchFlag = newFlag; //m_patchFlag is again just an unsigned integer
    	//auto just means that compiler assigns datatzype by itself
    	//const auto& mesh = solver.localMeshPtr(); // variable mesh which we use later for for loop; we assign a local Mesh pointer to it
    	const auto& mesh = solver.fullMeshPtr();


     	// Create patches by changing the markerID (flag) locally
    	unsigned int numNodesOnPatch(0); //here we just initalise an unsigned integer variable
    	for (int j(0); j < mesh->numBoundaryFacets(); j++) //returns number of boundary facets
    	{
    		auto& face = mesh->boundaryFacet(j);

                   for (int k(0); k < 3; ++k)
                    {
  	                      auto coord = face.point(k).coordinates(); //These cooridnates we want to get


  	                      auto currentprocessor = solver.comm()->MyPID();
  	                      std::ostringstream oss;//this is to convert int to string
  	                      oss << currentprocessor;

  	                      std::string path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinatefiles/coordinates_" + oss.str() + ".dat";
  	                      std::ofstream writer(path.c_str(), std::ios_base::app);
  	                      if(writer.is_open()==false)
  	                      {
  	                          	std::cout << "error occured while opening the file" << std::endl;
  	                      }
  	                      writer << coord[0] << "\t\t" << coord[1] << "\t\t" << coord[2] << std::endl;
  	                      writer.close();

   				          //auto pointInPatch = nodeOnPatch(coord, time); //here we check if the point is on the patch; if yes we set it to true

                    }

    	}


    }

*/




    vector_Type& patchDisplacementSum()
    {
        return *m_patchDisplacementVecSumPtr;
    }

    vectorPtr_Type patchDisplacementSumPtr()
    {
        return m_patchDisplacementVecSumPtr;
    }

    vector_Type& patchLocationSum()
    {
        return *m_patchLocationScalarSumPtr;
    }

    vectorPtr_Type patchLocationSumPtr()
    {
        return m_patchLocationScalarSumPtr;
    }

    vectorPtr_Type patchFacesLocationSumPtr()
    {
    	return m_patchFacesLocationScalarSumPtr;
    }

    vectorPtr_Type patchVecSumPtr()
    {
        	return m_patchVecSumPtr;
    }

    vectorPtr_Type directionalVecSumPtr()
    {
            	return m_directionVecFieldPtr;
    }


protected:
    
    void updatePatchDisplacementSum(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
        *m_patchDisplacementVecSumPtr *= 0.0;

        for (auto& patch : m_patchBCPtrVec)
        {
            *m_patchDisplacementVecSumPtr += patch->patchDisplacement(solver);
        }
    }
    
    void updatePatchLocationSum(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
        *m_patchLocationScalarSumPtr *= 0.0;

        for (auto& patch : m_patchBCPtrVec)
        {
            *m_patchLocationScalarSumPtr += patch->patchLocation();
        }
    }
    

    void updatePatchFacesLocationSum(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
    	*m_patchFacesLocationScalarSumPtr *= 0.0;
    	for (auto& patch : m_patchBCPtrVec)
    	{
    	            *m_patchFacesLocationScalarSumPtr += patch->patchFacesLocation();
    	}

    }

    void updatePatchVec(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {

    	*m_patchVecSumPtr *= 0.0;
    	for (auto& patch : m_patchBCPtrVec)
    	{
    	    	            *m_patchVecSumPtr += patch->patchVectorField(solver, time);
    	}
    }

    void updateDispDirectionalVec(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {
    	*m_directionVecFieldPtr *=0.0;

    	*m_directionVecFieldPtr += m_patchBCPtrVec[0]->displayDirectionalVectorField(solver, time);

    }


    const std::string m_patchListName;
    const GetPot& m_dataFile;
    const int m_patchNumber;

    vectorPtr_Type m_patchLocationScalarSumPtr;
    vectorPtr_Type m_patchDisplacementVecSumPtr;
    vectorPtr_Type m_patchFacesLocationScalarSumPtr;
    vectorPtr_Type m_patchVecSumPtr;
    vectorPtr_Type m_directionVecFieldPtr;

    std::vector<EssentialPatchBC*> m_patchBCPtrVec;

};
    
    
}

#endif /* EssentialPatchBC_hpp */

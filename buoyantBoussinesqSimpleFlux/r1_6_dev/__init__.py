#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV
##


#---------------------------------------------------------------------------   
def readTransportProperties( U, phi ):
    from Foam.transportModels import singlePhaseTransportModel
    laminarTransport = singlePhaseTransportModel( U, phi )
    
    from Foam.OpenFOAM import dimensionedScalar, word
    # Thermal expansion coefficient [1/K]
    beta = dimensionedScalar( laminarTransport.lookup( word( "beta" ) ) )

    # Reference temperature [K]
    TRef = dimensionedScalar( laminarTransport.lookup( word( "TRef" ) ) )

    # Laminar Prandtl number
    Pr = dimensionedScalar( laminarTransport.lookup( word( "Pr" ) ) )

    # Turbulent Prandtl number
    Prt = dimensionedScalar( laminarTransport.lookup( word( "Prt" ) ) )
    
    return laminarTransport, beta, TRef,Pr, Prt


#-------------------------------------------------------------------------
def createFields( runTime, mesh, g ):
    from Foam.OpenFOAM import ext_Info, nl
    from Foam.OpenFOAM import IOdictionary, IOobject, word, fileName
    from Foam.finiteVolume import volScalarField
    
    ext_Info() << "Reading thermophysical properties\n" << nl 

    ext_Info() << "Reading field T\n" << nl 
    T = volScalarField( IOobject( word( "T" ),
                                  fileName( runTime.timeName() ),
                                  mesh,
                                  IOobject.MUST_READ,
                                  IOobject.AUTO_WRITE ),
                        mesh )

    ext_Info() << "Reading field p\n" << nl
    p = volScalarField( IOobject( word( "p" ),
                                  fileName( runTime.timeName() ),
                                  mesh,
                                  IOobject.MUST_READ,
                                  IOobject.AUTO_WRITE ),
                        mesh )

    ext_Info() << "Reading field U\n" << nl
    from Foam.finiteVolume import volVectorField
    U = volVectorField( IOobject( word( "U" ),
                                  fileName( runTime.timeName() ),
                                  mesh,
                                  IOobject.MUST_READ,
                                  IOobject.AUTO_WRITE ),
                        mesh )
    
  
    from Foam.finiteVolume.cfdTools.incompressible import createPhi
    phi = createPhi( runTime, mesh, U )
    
    laminarTransport, beta, TRef,Pr, Prt = readTransportProperties( U, phi )
    
    ext_Info() << "Creating turbulence model\n" << nl
    from Foam import incompressible
    turbulence = incompressible.RASModel.New( U, phi, laminarTransport )
    
    ext_Info() << "Calculating field beta*(g.h)\n" << nl
    from Foam.finiteVolume import surfaceScalarField
    betaghf = surfaceScalarField( word( "betagh" ), beta * ( g & mesh.Cf() ) )
    
    pRefCell = 0
    pRefValue = 0.0
    
    from Foam.finiteVolume import setRefCell
    pRefCell, pRefValue = setRefCell( p, mesh.solutionDict().subDict( word( "SIMPLE" ) ), pRefCell, pRefValue )

    # Kinematic density for buoyancy force
    rhok = volScalarField( IOobject( word( "rhok" ),
                                     fileName( runTime.timeName() ),
                                     mesh ),
                           1.0 - beta * ( T - TRef ) )
    
    return T, p, U, phi, laminarTransport, beta, TRef,Pr, Prt, turbulence, betaghf, pRefCell, pRefValue, rhok


#--------------------------------------------------------------------------------------
def initConvergenceCheck( simple ):
    eqnResidual = 1
    maxResidual = 0
    convergenceCriterion = 0
    
    from Foam.OpenFOAM import word
    tmp, convergenceCriterion = simple.readIfPresent( word( "convergence" ), convergenceCriterion )
    
    return eqnResidual, maxResidual, convergenceCriterion


#--------------------------------------------------------------------------------------
def fun_UEqn( phi, U, p, turbulence, mesh, g, rhok, eqnResidual, maxResidual ):
    from Foam import fvm, fvc
    UEqn = fvm.div( phi, U ) + turbulence.divDevReff( U ) 
    UEqn.relax()
    
    from Foam.finiteVolume import solve
    eqnResidual = solve( UEqn == fvc.reconstruct( ( fvc.interpolate( rhok ) * ( g & mesh.Sf() ) - fvc.snGrad( p ) * mesh.magSf() ) ) ).initialResidual()
    
    maxResidual = max(eqnResidual, maxResidual)
       
    return UEqn, eqnResidual, maxResidual


#--------------------------------------------------------------------------------------
def fun_TEqn( turbulence, phi, T, rhok, beta, TRef, Pr, Prt, eqnResidual, maxResidual ):
    from Foam.OpenFOAM import word
    from Foam.finiteVolume import volScalarField
    kappaEff = volScalarField( word( "kappaEff" ),
                               turbulence.nu() / Pr + turbulence.ext_nut() / Prt )

    from Foam import fvc, fvm
    TEqn = fvm.div( phi, T ) - fvm.Sp( fvc.div( phi ), T ) - fvm.laplacian( kappaEff, T ) 

    TEqn.relax()

    eqnResidual = TEqn.solve().initialResidual()
    maxResidual = max(eqnResidual, maxResidual)

    rhok.ext_assign( 1.0 - beta * ( T - TRef ) )
    
    return TEqn, kappaEff
    

#--------------------------------------------------------------------------------------    
def fun_pEqn( runTime, mesh, p, phi, U, UEqn, g, rhok, eqnResidual, maxResidual, nNonOrthCorr, cumulativeContErr, pRefCell, pRefValue ): 
    
    from Foam.finiteVolume import volScalarField, surfaceScalarField
    from Foam.OpenFOAM import word
    from Foam import fvc
    rUA = volScalarField( word( "rUA" ), 1.0 / UEqn().A() )
    rUAf = surfaceScalarField(word( "(1|A(U))" ), fvc.interpolate( rUA ) )

    U.ext_assign( rUA * UEqn().H() )
    UEqn.clear()
    
    from Foam import fvc 
    phi.ext_assign( fvc.interpolate( U ) & mesh.Sf() )
    
    from Foam.finiteVolume import adjustPhi
    adjustPhi( phi, U, p )
    
    buoyancyPhi = rUAf * fvc.interpolate( rhok ) * ( g & mesh.Sf() )
    
    phi.ext_assign( phi + buoyancyPhi )

    for nonOrth in range( nNonOrthCorr+1 ):
        
        from Foam import fvm, fvc
        pEqn = fvm.laplacian(rUAf, p) == fvc.div(phi)

        pEqn.setReference( pRefCell, pRefValue )
        
        # retain the residual from the first iteration
        if ( nonOrth == 0 ):
              eqnResidual = pEqn.solve().initialResidual()
              maxResidual = max( eqnResidual, maxResidual )
              pass
        else:
              pEqn.solve()
              pass
        

        if ( nonOrth == nNonOrthCorr ):
           # Calculate the conservative fluxes
           phi.ext_assign( phi - pEqn.flux() )
           
           # Explicitly relax pressure for momentum corrector
           p.relax()

           # Correct the momentum source with the pressure gradient flux
           # calculated from the relaxed pressure
           U.ext_assign( U + rUA * fvc.reconstruct( ( buoyancyPhi - pEqn.flux() ) / rUAf ) )
           U.correctBoundaryConditions()
           pass
        
        pass

    from Foam.finiteVolume.cfdTools.incompressible import continuityErrs
    cumulativeContErr = continuityErrs( mesh, phi, runTime, cumulativeContErr )
    
    return eqnResidual, maxResidual, cumulativeContErr


#--------------------------------------------------------------------------------------
def convergenceCheck( maxResidual, convergenceCriterion ):
    from Foam.OpenFOAM import ext_Info, nl
    
    if (maxResidual < convergenceCriterion):
       ext_Info() << "reached convergence criterion: " << convergenceCriterion << nl
       runTime.writeAndEnd();
       ext_Info   << "latestTime = " << runTime.timeName() << nl
       pass

    
#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    from Foam.OpenFOAM.include import setRootCase
    args = setRootCase( argc, argv )

    from Foam.OpenFOAM.include import createTime
    runTime = createTime( args )

    from Foam.OpenFOAM.include import createMesh
    mesh = createMesh( runTime )

    from Foam.finiteVolume.cfdTools.general.include import readGravitationalAcceleration
    g = readGravitationalAcceleration( runTime, mesh)
    
    T, p, U, phi, laminarTransport, beta, TRef,Pr, Prt, turbulence, betaghf, pRefCell, pRefValue, rhok = createFields( runTime, mesh, g )
    
    from Foam.finiteVolume.cfdTools.general.include import initContinuityErrs
    cumulativeContErr = initContinuityErrs()

    from Foam.OpenFOAM import ext_Info, nl
    ext_Info() << "\nStarting time loop\n" <<nl
    
    while runTime.loop():
        ext_Info() << "Time = " << runTime.timeName() << nl << nl

        from Foam.finiteVolume.cfdTools.general.include import readSIMPLEControls
        simple, nNonOrthCorr, momentumPredictor, fluxGradp, transonic = readSIMPLEControls( mesh )
        
        eqnResidual, maxResidual, convergenceCriterion = initConvergenceCheck( simple )
                
        p.storePrevIter()
        
        UEqn, eqnResidual, maxResidual = fun_UEqn( phi, U, p, turbulence, mesh, g, rhok, eqnResidual, maxResidual )
        
        TEqn, kappaEff = fun_TEqn( turbulence, phi, T, rhok, beta, TRef, Pr, Prt, eqnResidual, maxResidual )
        eqnResidual, maxResidual, cumulativeContErr = fun_pEqn( runTime, mesh, p, phi, U, UEqn, g, rhok, eqnResidual, \
                                                            maxResidual, nNonOrthCorr, cumulativeContErr, pRefCell, pRefValue )
        
        turbulence.correct()

        runTime.write()
        
        ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << \
              "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << nl
        
        convergenceCheck( maxResidual, convergenceCriterion ) 
        
        pass
        
    ext_Info() << "End\n" << nl 

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
from Foam import FOAM_BRANCH_VERSION
import sys, os
if FOAM_BRANCH_VERSION( "dev", ">=", "010600" ):
   if __name__ == "__main__" :
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
else:
   from Foam.OpenFOAM import ext_Info
   ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam1.6-ext or higher \n "

    
#--------------------------------------------------------------------------------------


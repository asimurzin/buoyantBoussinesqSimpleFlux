#!/usr/bin/env python

#---------------------------------------------------------------------------
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
from Foam import man, ref


#---------------------------------------------------------------------------
def readGravitationalAcceleration( runTime, mesh ):
    ref.ext_Info() << "\nReading g" << ref.nl
    
    g= ref.uniformDimensionedVectorField( ref.IOobject( ref.word("g"),
                                                        ref.fileName( runTime.constant() ),
                                                        mesh,
                                                        ref.IOobject.MUST_READ,
                                                        ref.IOobject.NO_WRITE ) )
    return g


#---------------------------------------------------------------------------
def readTransportProperties( U, phi ):
  
  laminarTransport = man.singlePhaseTransportModel( U, phi )

  # Thermal expansion coefficient [1/K]
  beta = ref.dimensionedScalar( laminarTransport.lookup( ref.word( "beta" ) ) )

  # Reference temperature [K]
  TRef = ref.dimensionedScalar( laminarTransport.lookup( ref.word( "TRef" ) ) )

  # Laminar Prandtl number
  Pr = ref.dimensionedScalar( laminarTransport.lookup( ref.word( "Pr" ) ) )

  # Turbulent Prandtl number
  Prt = ref.dimensionedScalar( laminarTransport.lookup( ref.word( "Prt" ) ) )
  
  return laminarTransport, beta, TRef, Pr, Prt


#---------------------------------------------------------------------------
def createFields( runTime, mesh, g ):
  
  ref.ext_Info() << "Reading thermophysical properties\n" << ref.nl
  
  ref.ext_Info() << "Reading field T\n" << ref.nl
  T = man.volScalarField( man.IOobject( ref.word( "T" ),
                                        ref.fileName( runTime.timeName() ),
                                        mesh,
                                        ref.IOobject.MUST_READ,
                                        ref.IOobject.AUTO_WRITE ), mesh )

  ref.ext_Info() << "Reading field p_rgh\n" << ref.nl
  p_rgh = man.volScalarField( man.IOobject( ref.word( "p_rgh" ),
                                            ref.fileName( runTime.timeName() ),
                                            mesh,
                                            ref.IOobject.MUST_READ,
                                            ref.IOobject.AUTO_WRITE ),
                                mesh )
  
  ref.ext_Info() << "Reading field U\n" << ref.nl
  U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                        ref.fileName( runTime.timeName() ),
                                        mesh,
                                        ref.IOobject.MUST_READ,
                                        ref.IOobject.AUTO_WRITE ), mesh )
  phi = man.createPhi( runTime, mesh, U )
    
  laminarTransport, beta, TRef, Pr, Prt = readTransportProperties( U, phi )
  
  ref.ext_Info()<< "Creating turbulence model\n" << ref.nl
  turbulence = man.incompressible.RASModel.New(U, phi, laminarTransport)

  # Kinematic density for buoyancy force
  rhok = man.volScalarField( man.IOobject( ref.word( "rhok" ),
                                           ref.fileName( runTime.timeName() ),
                                           mesh ), man( 1.0 - beta * ( T() - TRef ), man.Deps( T ) ) )
  
  # kinematic turbulent thermal thermal conductivity m2/s
  ref.ext_Info() << "Reading field kappat\n" << ref.nl
  kappat = man.volScalarField( man.IOobject( ref.word( "kappat" ),
                                             ref.fileName( runTime.timeName() ),
                                             mesh,
                                             ref.IOobject.MUST_READ,
                                             ref.IOobject.AUTO_WRITE ), mesh )

  ref.ext_Info() << "Calculating field g.h\n" << ref.nl
  gh = man.volScalarField( ref.word( "gh" ), man( g & mesh.C(), man.Deps( mesh ) ) )
  
  ghf = man.surfaceScalarField( ref.word( "ghf" ), man( g & mesh.Cf(), man.Deps( mesh ) ) )

  p = man.volScalarField( man.IOobject( ref.word( "p" ),
                                        ref.fileName( runTime.timeName() ),
                                        mesh,
                                        ref.IOobject.NO_READ,
                                        ref.IOobject.AUTO_WRITE ), p_rgh + rhok * gh )

  pRefCell = 0
  pRefValue = 0.0

  pRefCell, pRefValue = ref.setRefCell( p, p_rgh, mesh.solutionDict().subDict( ref.word( "SIMPLE" ) ), pRefCell, pRefValue )

  if p_rgh.needReference():
    p().ext_assign( p() + ref.dimensionedScalar( ref.word( "p" ),p.dimensions(), pRefValue - ref.getRefCellValue( p, pRefCell ) ) )
    pass
  
  return T, p_rgh, U, phi, laminarTransport, turbulence, rhok, kappat, gh, ghf, p, pRefCell, pRefValue, beta, TRef, Pr, Prt


#---------------------------------------------------------------------------
def fun_UEqn( mesh, simple, U, phi, turbulence, p, rhok, p_rgh, ghf ):
   
  UEqn = man.fvm.div(phi, U) + man( turbulence.divDevReff( U ) , man.Deps( U ) ) 

  UEqn.relax()

  if  simple.momentumPredictor():
    ref.solve( UEqn == man.fvc.reconstruct( ( - ghf * man.fvc.snGrad( rhok ) - man.fvc.snGrad( p_rgh ) ) * man.surfaceScalarField( mesh.magSf(), mesh ) ) )
    pass
  return UEqn


#---------------------------------------------------------------------------
def fun_TEqn( phi, turbulence, kappat, T, rhok, beta, TRef, Prt, Pr ):
  
  kappat().ext_assign( turbulence.ext_nut() / Prt )
  
  kappat.correctBoundaryConditions()
  
  kappaEff = ref.volScalarField ( ref.word( "kappaEff" ) , turbulence.ext_nu() / Pr  + kappat )
  TEqn = ref.fvm.div( phi, T ) - ref.fvm.Sp( ref.fvc.div( phi ), T ) - ref.fvm.laplacian( kappaEff, T ) 

  TEqn.relax()
  TEqn.solve()

  rhok().ext_assign( 1.0 - beta * ( T() - TRef ) )
  pass


#---------------------------------------------------------------------------
def fun_pEqn( mesh, runTime, simple, p, rhok, U, phi, turbulence, gh, ghf, p_rgh, UEqn, pRefCell, pRefValue, cumulativeContErr ):

  rAU = ref.volScalarField( ref.word( "rAU" ), 1.0 / UEqn.A() )
  
  rAUf = ref.surfaceScalarField( ref.word( "(1|A(U))" ), ref.fvc.interpolate( rAU ) )

  U().ext_assign( rAU * UEqn.H() )

  phi().ext_assign( ref.fvc.interpolate( U ) & mesh.Sf() )

  ref.adjustPhi( phi, U, p_rgh );

  buoyancyPhi = rAUf * ghf() * ref.fvc.snGrad( rhok ) * mesh.magSf()
  phi().ext_assign( phi() - buoyancyPhi )

  for nonOrth in range( simple.nNonOrthCorr() + 1 ):
    p_rghEqn = ref.fvm.laplacian( rAUf, p_rgh ) == ref.fvc.div( phi )

    p_rghEqn.setReference( pRefCell, ref.getRefCellValue( p_rgh, pRefCell ) )
    
    p_rghEqn.solve()

    if nonOrth == simple.nNonOrthCorr():
      # Calculate the conservative fluxes
      phi().ext_assign( phi() - p_rghEqn.flux() )
      
      # Explicitly relax pressure for momentum corrector
      p_rgh.relax()

      # Correct the momentum source with the pressure gradient flux
      # calculated from the relaxed pressure
      U().ext_assign( U() - rAU * ref.fvc.reconstruct( ( buoyancyPhi + p_rghEqn.flux() ) / rAUf ) )
      U.correctBoundaryConditions()
      pass
    pass

  cumulativeContErr = ref.ContinuityErrs( phi, runTime, mesh, cumulativeContErr )

  p.ext_assign( p_rgh + rhok * gh )

  if p_rgh.needReference():
    p().ext_assign( p() + ref.dimensionedScalar( ref.word( "p" ), p.dimensions(), pRefValue - ref.getRefCellValue( p, pRefCell ) ) )
    p_rgh.ext_assign( p - rhok * gh )
    pass
  
  return cumulativeContErr



#---------------------------------------------------------------------------
def main_standalone( argc, argv ):
  
  args = ref.setRootCase( argc, argv )
  
  runTime=man.createTime( args )
    
  mesh = man.createMesh( runTime )
    
  g = readGravitationalAcceleration( runTime, mesh );

  T, p_rgh, U, phi, laminarTransport, turbulence, rhok, \
     kappat, gh, ghf, p, pRefCell, pRefValue, beta, TRef, Pr, Prt  = createFields( runTime, mesh, g )
                                                        
  cumulativeContErr = ref.initContinuityErrs()

  simple = man.simpleControl( mesh )


  # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

  ref.ext_Info() << "\nStarting time loop\n" << ref.nl

  while simple.loop():
   
    ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl

    p_rgh.storePrevIter()

    # Pressure-velocity SIMPLE corrector
    UEqn = fun_UEqn( mesh, simple, U, phi, turbulence, p, rhok, p_rgh, ghf )
    fun_TEqn( phi, turbulence, kappat, T, rhok, beta, TRef, Prt, Pr )

    cumulativeContErr = fun_pEqn( mesh, runTime, simple, p, rhok, U, phi, turbulence, gh, ghf, p_rgh, UEqn, pRefCell, pRefValue, cumulativeContErr )
    
    turbulence.correct()
    runTime.write()
    
    ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" \
               << "  ClockTime = " << runTime.elapsedClockTime() << " s" \
               << ref.nl << ref.nl
    pass
  pass

  ref.ext_Info() << "End\n" << ref.nl

  import os
  return os.EX_OK


#---------------------------------------------------------------------------
from Foam import FOAM_REF_VERSION
if FOAM_REF_VERSION( ">=", "020000" ):
   if __name__ == "__main__" :
      import sys, os
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
   ref.ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam2.0.0 \n "     
   pass


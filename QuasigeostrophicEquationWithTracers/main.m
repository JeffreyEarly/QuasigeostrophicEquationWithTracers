//
//  main.m
//  QuasigeostrophicEquationWithTracers
//
//  Created by Jeffrey Early on 4/13/12.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

//#import <GLNumericalModelingKit/GLOperationOptimizer.h>

int main (int argc, const char * argv[])
{
	
	@autoreleasepool {
	    
		// Reasonable parameters to nondimensionalize by.
		GLFloat N_QG = 1.3; // cm
		GLFloat T_QG = 12; // days
		GLFloat L_QG = 47; // km
		
		// 256x128 at takes 11 second with optimized code.
		GLDimension *xDim = [[GLDimension alloc] initPeriodicDimension: YES nPoints: 512 domainMin: -1500/L_QG length: 2000/L_QG];
		xDim.name = @"x";
		GLDimension *yDim = [[GLDimension alloc] initPeriodicDimension: YES nPoints: 256 domainMin: -500/L_QG length: 1000/L_QG];
		yDim.name = @"y";
		GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: [NSArray arrayWithObject: [NSNumber numberWithDouble: 0.0]]];
		tDim.name = @"time";
		
		NSLog(@"Spatial resolution: %f km", xDim.sampleInterval*L_QG);
		
		// Variables are always tied to a particular equation---so we create an equation object first.
		GLEquation *equation = [[GLEquation alloc] init];
		
		NSArray *spatialDimensions = [NSArray arrayWithObjects: xDim, yDim, nil];
		GLVariable *x = [GLVariable variableOfRealTypeFromDimension: xDim withDimensions: spatialDimensions forEquation: equation];
		GLVariable *y = [GLVariable variableOfRealTypeFromDimension: yDim withDimensions: spatialDimensions forEquation: equation];
		
		/************************************************************************************************/
		/*		Create and cache the differential operators that we will be using	                 	*/
		/************************************************************************************************/
		
		// At the moment we know that this is the spectral operators, although in the future we'll have to set this up explicitly.
		GLSpectralDifferentialOperatorPool *diffOperators = [equation defaultDifferentialOperatorPoolForVariable: x];
		
		// Create the operator xx+yy-1---this is how you compute y from eta
		GLSpectralDifferentialOperator *laplacianMinusOne = [[diffOperators harmonicOperator] scalarAdd: -1.0];
		[diffOperators setDifferentialOperator: laplacianMinusOne forName: @"laplacianMinusOne"];
		
		// Create the operator 1/(xx+yy-1)---this is how you compute eta from y.
		[diffOperators setDifferentialOperator: [laplacianMinusOne scalarDivide: 1.0] forName: @"inverseLaplacianMinusOne"];
		
		// This builds the differentiation matrix diff_{xxx} + diff_{xyy}
		[diffOperators setDifferentialOperator: [[diffOperators xxx] plus: [diffOperators xyy]] forName: @"diffJacobianX"];
		
		// This builds the differentiation matrix diff_{xxy} + diff_{yyy}
		[diffOperators setDifferentialOperator: [[diffOperators xxy] plus: [diffOperators yyy]] forName: @"diffJacobianY"];
		
		GLSpectralDifferentialOperator *svv = [diffOperators spectralVanishingViscosityFilter];
		
		GLFloat k = 0.05*xDim.sampleInterval;
		//GLSpectralDifferentialOperator *diffLin = [[[[diffOperators harmonicOperatorOfOrder: 2] scalarMultiply: k] multiply: svv] minus: [diffOperators x]];
		GLSpectralDifferentialOperator *diffLin = [[[diffOperators harmonicOperatorOfOrder: 2] scalarMultiply: k] minus: [diffOperators x]];
		[diffOperators setDifferentialOperator: diffLin forName: @"diffLin"];
		
		GLSpectralDifferentialOperator *harmonicDamp = [[[diffOperators harmonicOperatorOfOrder: 1] scalarMultiply: k] multiply: svv];
		[diffOperators setDifferentialOperator: harmonicDamp forName: @"damp"];
		
		/************************************************************************************************/
		/*		Create the initial conditions for the ssh, tracer, and floats							*/
		/************************************************************************************************/
		
		//  gaussian = amplitude * exp( - ((x-x0)*(x-x0) + (y-y0)*(y-y0))/(length*length) );
		GLFloat amplitude = 15.0/N_QG;
		GLFloat length = 80/L_QG;
		
		GLVariable *r2 = [[x times: x] plus: [y times: y]];
		GLVariable *gaussian = [[[r2 scalarMultiply: -1.0/(length*length)] exponentiate] scalarMultiply: amplitude];
		
		// The tracer varies linearly from east-to-west, but we create a smooth edge at the boundary.
		GLVariable *distanceFromEW = [[[[[x scalarAdd: -xDim.domainMin-xDim.domainLength/2.0] abs] negate] scalarAdd: xDim.domainLength/2.0] scalarMultiply: 1/(250/L_QG)];
		GLVariable *tracer = [[[[[[distanceFromEW times: distanceFromEW] scalarMultiply: -2.0*M_PI] exponentiate] negate] scalarAdd: 1] times: x];
		
		// Let's also plop a float at each grid point.
		GLVariable *xPosition = [GLVariable variableFromVariable: x];
		GLVariable *yPosition = [GLVariable variableFromVariable: y];
        
		//        GLDimension *drifterDim = [[GLDimension alloc] initPeriodicDimension: NO nPoints:xDim.nPoints domainMin:1 sampleInterval:1];
		//		drifterDim.name = @"drifterID";
		//		GLVariable *xPosition = [GLVariable variableOfRealTypeWithDimensions: [NSArray arrayWithObject: drifterDim] forEquation:equation];
		//		GLVariable *yPosition = [GLVariable variableOfRealTypeWithDimensions: [NSArray arrayWithObject: drifterDim] forEquation:equation];
		//        for ( NSUInteger i=0; i<drifterDim.nPoints; i++){
		//            xPosition.pointerValue[i] = [xDim valueAtIndex: i];
		//            yPosition.pointerValue[i] = 0;
		//        }
		
		//        xPosition.pointerValue[0] = [xDim valueAtIndex: 364];
		//        yPosition.pointerValue[1] = 0;
		//
		//        xPosition.pointerValue[0] = [xDim valueAtIndex: 424];
		//        yPosition.pointerValue[1] = 0;
        
		/************************************************************************************************/
		/*		Create a NetCDF file and mutable variables in order to record some of the time steps.	*/
		/************************************************************************************************/
		
		GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: @"/Users/jearly/Desktop/QuasigeostrophyTracers.nc"] forEquation: equation overwriteExisting: YES];
		[netcdfFile setGlobalAttribute: @(N_QG) forKey: @"height_scale"];
		[netcdfFile setGlobalAttribute: @(L_QG) forKey: @"length_scale"];
		[netcdfFile setGlobalAttribute: @(T_QG) forKey: @"time_scale"];
		GLMutableVariable *sshHistory = [gaussian variableByAddingDimension: tDim];
		sshHistory.name = @"SSH";
		sshHistory = [netcdfFile addVariable: sshHistory];
		GLMutableVariable *tracerHistory = [tracer variableByAddingDimension: tDim];
		tracerHistory.name = @"x-tracer";
		tracerHistory = [netcdfFile addVariable: tracerHistory];
		GLMutableVariable *xPositionHistory = [xPosition variableByAddingDimension: tDim];
		xPositionHistory.name = @"x-position";
		xPositionHistory = [netcdfFile addVariable: xPositionHistory];
		GLMutableVariable *yPositionHistory = [yPosition variableByAddingDimension: tDim];
		yPositionHistory.name = @"y-position";
		yPositionHistory = [netcdfFile addVariable: yPositionHistory];
		
		/************************************************************************************************/
		/*		Determine an appropriate time step based on the CFL condition.							*/
		/************************************************************************************************/
		
		GLVariable *v = [gaussian x];
		GLVariable *u = [gaussian y];
		GLVariable *speed = [[u times: u] plus: [v times: v]];
		[equation solveForVariable: speed];
		
		CGFloat cfl = 0.5;
		GLFloat U = sqrt([speed maxNow]);
		GLFloat timeStep = cfl * xDim.sampleInterval / U;
		GLFloat maxTime = 365/T_QG;
		
		/************************************************************************************************/
		/*		Create the integration object.															*/
		/************************************************************************************************/
		
		GLVariable *ssh = [gaussian diff: @"laplacianMinusOne"];
		tracer = [tracer frequencyDomain];
		NSArray *yin = @[ssh, tracer, xPosition, yPosition];
		
		GLVectorIntegrationOperation *integrator = [GLVectorIntegrationOperation rungeKutta4AdvanceY: yin stepSize: timeStep fFromY:^(NSArray *yNew) {
			//		GLAdaptiveVectorIntegrationOperation *integrator = [GLAdaptiveVectorIntegrationOperation rungeKutta45AdvanceY: yin stepSize: timeStep fFromY:^(NSArray *yNew) {
			GLVariable *eta = [yNew[0] diff: @"inverseLaplacianMinusOne"];
			
			GLVariable *fSSH = [[eta diff:@"diffLin"] plus: [[[[eta y] times: [eta diff: @"diffJacobianX"]] minus: [[eta x] times: [eta diff: @"diffJacobianY"]]] frequencyDomain]];
			
			GLVariable *yTracer = yNew[1];
			GLVariable *fTracer = [[[[[eta y] times: [yTracer x]] minus:[[eta x] times: [yTracer y]] ] frequencyDomain] plus: [yTracer diff: @"damp"]];
			
			NSArray *uv = @[[[[eta y] spatialDomain] negate], [[eta x] spatialDomain] ];
			NSArray *xy = @[yNew[2], yNew[3]];
			GLInterpolationOperation *interp = [[GLInterpolationOperation alloc] initWithFirstOperand: uv secondOperand: xy];
			
			NSArray *f = @[fSSH, fTracer, interp.result[0], interp.result[1]];
			return f;
		}];
		
		//		integrator.absoluteTolerance = @[ @(1e-6), @(1e-6), @(1e-3), @(1e-3)];
		
		/************************************************************************************************/
		/*		Step forward in time, and write the data to file every-so-often.						*/
		/************************************************************************************************/
		
		for (GLFloat time = 1/T_QG; time < maxTime; time += 1/T_QG)
		{
            @autoreleasepool {
				yin = [integrator stepForward: yin toTime: time];
				
				ssh = yin[0];
				tracer = yin[1];
				xPosition = yin[2];
				yPosition = yin[3];
				
				NSLog(@"Logging day: %f, step size: %f.", (integrator.currentTime*T_QG), integrator.lastStepSize*T_QG);
				// We're using spectral code, so it's possible (and is in fact the case) that the variable is not in the spatial domain.
				[tDim addPoint: [NSNumber numberWithDouble: time]];
				GLVariable *eta = [[ssh diff: @"inverseLaplacianMinusOne"] spatialDomain];
				[sshHistory concatenateWithLowerDimensionalVariable: eta alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[tracerHistory concatenateWithLowerDimensionalVariable: [tracer spatialDomain] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[xPositionHistory concatenateWithLowerDimensionalVariable: xPosition alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[yPositionHistory concatenateWithLowerDimensionalVariable: yPosition alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
            }
		}
		
		[equation waitUntilAllOperationsAreFinished];
		
		// The NetCDF file may still be writing data. We need to make sure it finishes before we exit.
		[netcdfFile waitUntilAllOperationsAreFinished];
		[netcdfFile close];
		
	    
	}
    return 0;
}


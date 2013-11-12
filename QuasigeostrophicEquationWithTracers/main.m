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
		
		/************************************************************************************************/
		/*		Define the problem dimensions															*/
		/************************************************************************************************/
		
		GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:512 domainMin:-1500/L_QG length:2000/L_QG];
		xDim.name = @"x";
		GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:256 domainMin:-500/L_QG length:1000/L_QG];
		yDim.name = @"y";
		GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
		tDim.name = @"time";
		
		// Variables are always tied to a particular equation---so we create an equation object first.
		GLEquation *equation = [[GLEquation alloc] init];
		
		NSArray *spatialDimensions = @[xDim, yDim];
		GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: spatialDimensions forEquation: equation];
		GLFunction *y = [GLFunction functionOfRealTypeFromDimension: yDim withDimensions: spatialDimensions forEquation: equation];
		
		/************************************************************************************************/
		/*		Create and cache the differential operators we will need								*/
		/************************************************************************************************/
		
		NSArray *spectralDimensions = [x dimensionsTransformedToBasis: x.differentiationBasis];
		
		GLLinearTransform *laplacian = [GLLinearTransform harmonicOperatorFromDimensions: spectralDimensions forEquation: equation];
		GLLinearTransform *laplacianMinusOne = [laplacian plus: @(-1.0)];
		GLLinearTransform *inverseLaplacianMinusOne = [laplacianMinusOne inverse];
		
		GLLinearTransform *diff_xxx = [GLLinearTransform differentialOperatorWithDerivatives:@[@(3),@(0)] fromDimensions:spectralDimensions forEquation:equation];
		GLLinearTransform *diff_xyy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(1),@(2)] fromDimensions:spectralDimensions forEquation:equation];
		GLLinearTransform *diff_xxy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(2),@(1)] fromDimensions:spectralDimensions forEquation:equation];
		GLLinearTransform *diff_yyy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(0),@(3)] fromDimensions:spectralDimensions forEquation:equation];
		
		GLLinearTransform *diffJacobianX = [diff_xxx plus: diff_xyy];
		GLLinearTransform *diffJacobianY = [diff_xxy plus: diff_yyy];
		
		GLFloat k = 0.05*xDim.sampleInterval;
		GLLinearTransform *biharmonic = [GLLinearTransform harmonicOperatorOfOrder: 2 fromDimensions: spectralDimensions forEquation: equation];
		GLLinearTransform *svv = [GLLinearTransform spectralVanishingViscosityFilterWithDimensions: spectralDimensions scaledForAntialiasing: YES forEquation: equation];
		GLLinearTransform *diffLin = [[biharmonic times: @(k)] times: svv];
		
		GLLinearTransform *harmonicDamp = [[laplacian times: @(k)] times: svv];
		
		/************************************************************************************************/
		/*		Create the initial conditions for the ssh, tracer, and floats							*/
		/************************************************************************************************/
		
		//  gaussian = amplitude * exp( - ((x-x0)*(x-x0) + (y-y0)*(y-y0))/(length*length) );
		GLFloat amplitude = 15.0/N_QG;
		GLFloat length = 80/L_QG;
		
		GLFunction *r2 = [[x times: x] plus: [y times: y]];
		GLFunction *gaussian = [[[r2 times: @(-1.0/(length*length))] exponentiate] times: @(amplitude)];
		
		// The tracer varies linearly from east-to-west, but we create a smooth edge at the boundary.
		GLFunction *distanceFromEW = [[[[[x plus: @(-xDim.domainMin-xDim.domainLength/2.0)] abs] negate] plus: @(xDim.domainLength/2.0)] times: @(1/(250/L_QG))];
		GLFunction *tracer = [[[[[[distanceFromEW times: distanceFromEW] times: @(-2.0*M_PI)] exponentiate] negate] plus: @(1)] times: x];
		
		// Let's also plop a float at each grid point.
		GLFunction *xPosition = [GLFunction functionFromFunction: x];
		GLFunction *yPosition = [GLFunction functionFromFunction: y];
        
		/************************************************************************************************/
		/*		Create a NetCDF file and mutable variables in order to record some of the time steps.	*/
		/************************************************************************************************/
		
		NSString *path = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent:@"QuasigeostrophyTracers.nc"];
		GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: path] forEquation: equation overwriteExisting: YES];
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
		
		GLFunction *v = [gaussian x];
		GLFunction *u = [gaussian y];
		GLFunction *speed = [[u times: u] plus: [v times: v]];
		[equation solveForVariable: speed];
		
		CGFloat cfl = 0.5;
		GLFloat U = sqrt([speed maxNow]);
		GLFloat timeStep = cfl * xDim.sampleInterval / U;
		GLFloat maxTime = 365/T_QG;
		
		/************************************************************************************************/
		/*		Create the integration object.															*/
		/************************************************************************************************/
		
		GLFunction *ssh = [gaussian differentiateWithOperator: laplacianMinusOne];
		tracer = [tracer frequencyDomain];
		NSArray *yin = @[ssh, tracer, xPosition, yPosition];
		
//		GLRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta23AdvanceY: yin stepSize: timeStep fFromTY:^(GLScalar *time, NSArray *yNew) {
		GLRungeKuttaOperation *integrator = [GLRungeKuttaOperation rungeKutta4AdvanceY: yin stepSize: timeStep fFromTY:^(GLScalar *time, NSArray *yNew) {
			GLFunction *eta = [inverseLaplacianMinusOne transform: yNew[0]];
			
			GLFunction *fSSH = [[eta differentiateWithOperator: diffLin] plus: [[[[eta y] times: [eta differentiateWithOperator: diffJacobianX]] minus: [[eta x] times: [[[eta differentiateWithOperator: diffJacobianY] spatialDomain] plus: @(1.0)]]] frequencyDomain]];
			
			GLFunction *yTracer = yNew[1];
			GLFunction *fTracer = [[[[[eta y] times: [yTracer x]] minus:[[eta x] times: [yTracer y]] ] frequencyDomain] plus: [yTracer differentiateWithOperator: harmonicDamp]];
			
			NSArray *uv = @[[[[eta y] spatialDomain] negate], [[eta x] spatialDomain] ];
			NSArray *xy = @[yNew[2], yNew[3]];
			GLSimpleInterpolationOperation *interp = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: uv secondOperand: xy];
			
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
				yin = [integrator stepForwardToTime: time];
				
				ssh = yin[0];
				tracer = yin[1];
				xPosition = yin[2];
				yPosition = yin[3];
				
				NSLog(@"Logging day: %f, step size: %f.", (integrator.currentTime*T_QG), integrator.lastStepSize*T_QG);
				// We're using spectral code, so it's possible (and is in fact the case) that the variable is not in the spatial domain.
				[tDim addPoint: @(time)];
				GLFunction *eta = [[inverseLaplacianMinusOne transform: ssh] spatialDomain];
				GLFunction *tracerSpace = [tracer spatialDomain];
				[sshHistory concatenateWithLowerDimensionalVariable: eta alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[tracerHistory concatenateWithLowerDimensionalVariable: tracerSpace alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
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


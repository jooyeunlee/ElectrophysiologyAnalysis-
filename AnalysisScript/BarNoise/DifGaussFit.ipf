#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Written by David B. Kastner

function getDifGaussFit(wv, multiplier)
	variable multiplier
	wave wv
	K0 = 0;
	duplicate /o wv posVals
	posVals=limit(posVals,0,inf)
	CurveFit/q/w=2/H="1000"/NTHR=0 gauss posVals
	wave w_coef 
	Make/D/N=5/O params
	params[0] = {0,w_coef[3]/sqrt(2),w_coef[3]/sqrt(2)*multiplier,w_coef[1],w_coef[1]/2}
	variable peakPnt=(w_coef[2]-leftx(wv))/deltax(wv)

	make /o/n=(min(floor(peakPnt),numpnts(posVals)-ceil(peakPnt)-1)) halfSpace
	setscale /p x,0,deltax(wv),halfSpace
	halfSpace=wv(x+w_coef[2])
	halfSpace+=wv(w_coef[2]-x)
	
	Make/O/T/N=8 T_Constraints
	T_Constraints[0] = {"K1 > 1e-6","K1 < 20","K2 > 1e-6","K2 < 100","K3 > 0.05","K3 < 10","K4 > 0.0001","K4 < 10"}
	variable v_fitError=0
	FuncFit/q/NTHR=0/H="10000" difGauss params  halfSpace /D /C=T_Constraints
	variable noFit=0
	if(v_fitError>0)
		K0 = 0;K2=0;
		CurveFit/q/w=2/H="1010"/NTHR=0 gauss halfSpace /D
		duplicate /o w_coef params
		noFit=1
	endif
	
	killWindow /z gaussFitGraph
	doGaussFitGraph(noFit)
end

//for spike
function getDifGaussFitsp(wv, multiplier)
	variable multiplier
	wave wv
	K0 = 0;
	duplicate /o wv posVals
	posVals=limit(posVals,0,inf)
	CurveFit/q/w=2/H="1000"/NTHR=0 gauss posVals
	wave w_coef 
	Make/D/N=5/O params
	params[0] = {0,w_coef[3]/sqrt(2),w_coef[3]/sqrt(2)*multiplier,w_coef[1],w_coef[1]/2}
	variable peakPnt=(w_coef[2]-leftx(wv))/deltax(wv)

	make /o/n=(min(floor(peakPnt),numpnts(posVals)-ceil(peakPnt)-1)) halfSpace
	setscale /p x,0,deltax(wv),halfSpace
	halfSpace=wv(x+w_coef[2])
	halfSpace+=wv(w_coef[2]-x)
	
	Make/O/T/N=8 T_Constraints
	T_Constraints[0] = {"K1 > 1e-6","K1 < 20","K2 > 1e-6","K2 < 100","K3 > 0.05","K3 < 10","K4 > 0.0001","K4 < 10"}
	variable v_fitError=0
	FuncFit/q/NTHR=0/H="10000" difGauss params  halfSpace /D /C=T_Constraints
	variable noFit=0
	if(v_fitError>0)
		K0 = 0;K2=0;
		CurveFit/q/w=2/H="1010"/NTHR=0 gauss halfSpace /D
		duplicate /o w_coef params
		noFit=1
	endif
	
	killWindow /z gaussFitGraph
	doGaussFitGraph(noFit)
end

//for membrane voltage
function getDifGaussFitmem(wv, multiplier)
	variable multiplier
	wave wv
	K0 = 0;
	duplicate /o wv posVals
	posVals=limit(posVals,0,inf)
	CurveFit/q/w=2/H="1000"/NTHR=0 gauss posVals
	wave w_coef 
	Make/D/N=5/O paramsmem
	paramsmem[0] = {0,w_coef[3]/sqrt(2),w_coef[3]/sqrt(2)*multiplier,w_coef[1],w_coef[1]/2}
	variable peakPnt=(w_coef[2]-leftx(wv))/deltax(wv)

	make /o/n=(min(floor(peakPnt),numpnts(posVals)-ceil(peakPnt)-1)) halfSpace
	setscale /p x,0,deltax(wv),halfSpace
	halfSpace=wv(x+w_coef[2])
	halfSpace+=wv(w_coef[2]-x)
	
	Make/O/T/N=8 T_Constraints
	T_Constraints[0] = {"K1 > 1e-6","K1 < 20","K2 > 1e-6","K2 < 100","K3 > 0.05","K3 < 10","K4 > 0.0001","K4 < 10"}
	variable v_fitError=0
	FuncFit/q/NTHR=0/H="10000" difGauss paramsmem  halfSpace /D /C=T_Constraints
	variable noFit=0
	if(v_fitError>0)
		K0 = 0;K2=0;
		CurveFit/q/w=2/H="1010"/NTHR=0 gauss halfSpace /D
		duplicate /o w_coef paramsmem
		noFit=1
	endif
	
	killWindow /z gaussFitGraph
	doGaussFitGraph(noFit)
end

Function difGauss(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = A1*exp(-(x-mn)^2/2/sd1^2)- A2*exp(-(x-mn)^2/2/sd2^2)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = mn
	//CurveFitDialog/ w[1] = sd1
	//CurveFitDialog/ w[2] = sd2
	//CurveFitDialog/ w[3] = A1
	//CurveFitDialog/ w[4] = A2

	return w[3]*exp(-(x-w[0])^2/2/w[1]^2)- w[4]*exp(-(x-w[0])^2/2/w[2]^2)
End

function doGaussFitGraph(which)
	variable which
	wave halfSpace,fit_halfSpace
	Display /k=1/n=gaussFitGraph/W=(773,567,930,697) halfSpace,fit_halfSpace
	ModifyGraph lSize=2
	ModifyGraph rgb(halfSpace)=(0,0,0)
	if(which)
		ModifyGraph rgb(fit_halfSpace)=(65535,0,0)
	else
		ModifyGraph rgb(fit_halfSpace)=(3,52428,1)
	endif
	ModifyGraph zero(left)=2
	ModifyGraph nticks=2
	ModifyGraph fSize=10
	ModifyGraph standoff=0
	ModifyGraph ZisZ=1
	ModifyGraph zapTZ=1
	ModifyGraph btLen=1.5
	Label left "\\Z12\\F'Helvetica'Magnitude"
	Label bottom "\\Z12\\F'Helvetica'Size (mm)"
end

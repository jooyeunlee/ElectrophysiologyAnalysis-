#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Written by David B Kastner
// Edited by Joo Yeun Lee

CONSTANT NumBootstrap=10000

function getDimensions(cWv,hWv,dWv)
	wave cWv,hWv,dWv
	make /o/n=(3,dimsize(cWv,1)) pValSet
	make /o/n=(3,dimsize(cWv,1)) dimVals
	make /o/n=(3,dimsize(cWv,1)) stdevs 
	variable i,j
	for(i=0;i<dimsize(pValSet,1);i+=1)
		make /o/n=(dimsize(cWv,0)) c1=cWv[p][i]
		make /o/n=(dimsize(hWv,0)) h1=hWv[p][i]
		make /o/n=(dimsize(dWv,0)) d1=dWv[p][i]
		wavetransform zapNaNs,c1
		wavetransform zapNaNs,h1
		wavetransform zapNaNs,d1
	
		dimVals[0][i]=getMag(c1,h1,d1)
		dimVals[1][i]=getMag(h1,c1,d1)
		dimVals[2][i]=getMag(d1,c1,h1)
		
		wavestats /q c1
		stdevs[0][i]=V_sdev
		wavestats /q h1
		stdevs[1][i]=V_sdev
		wavestats /q d1
		stdevs[2][i]=V_sdev		
		
		statswilcoxonranktest /q/tail=4 c1,h1
		wave w_wilcoxonTest
		pValSet[0][i]=w_wilcoxonTest[5]
		statswilcoxonranktest /q/tail=4 c1,d1
		pValSet[1][i]=w_wilcoxonTest[5]
		statswilcoxonranktest /q/tail=4 h1,d1
		pValSet[2][i]=w_wilcoxonTest[5]
	endfor
	
	make /o/n=(dimsize(cWv,0)) c1=cWv[p][0],c2=cWv[p][1]
	make /o/n=(dimsize(hWv,0)) h1=hWv[p][0],h2=hWv[p][1]
	make /o/n=(dimsize(dWv,0)) d1=dWv[p][0],d2=dWv[p][1]
	wavetransform zapNaNs,c1
	wavetransform zapNaNs,h1
	wavetransform zapNaNs,d1
	wavetransform zapNaNs,c2
	wavetransform zapNaNs,h2
	wavetransform zapNaNs,d2
	duplicate /o c1 allC,holdC1
	duplicate /o c2 holdC2
	concatenate /NP=0 "c2;",allC
	duplicate /o h1 allH,holdH1
	duplicate /o h2 holdH2
	concatenate /NP=0 "h2;",allH
	duplicate /o d1 allD,holdD1
	duplicate /o d2 holdD2
	concatenate /NP=0 "d2;",allD
	
	make /o/n=(NumBootstrap,3) bootVals
	make /o/n=(3) pValAll,valAll
	for(i=0;i<NumBootstrap;i+=1)
		doBootAll(holdC1,holdC2,allC)
		doBootAll(holdH1,holdH2,allH)
		doBootAll(holdD1,holdD2,allD)
		bootVals[i][0]=abs(getMag(holdC1,holdH1,holdD1)-getMag(holdC2,holdH2,holdD2))
		bootVals[i][1]=abs(getMag(holdH1,holdC1,holdD1)-getMag(holdH2,holdC2,holdD2))
		bootVals[i][2]=abs(getMag(holdD1,holdC1,holdH1)-getMag(holdD2,holdC2,holdH2))
	endfor

	valAll[0]=getMag(c1,h1,d1)-getMag(c2,h2,d2)
	valAll[1]=getMag(h1,c1,d1)-getMag(h2,c2,d2)
	valAll[2]=getMag(d1,c1,h1)-getMag(d2,c2,h2)
	for(j=0;j<dimsize(pValAll,0);j+=1)
		make /o/n=(dimsize(bootVals,0)) hold=bootVals[p][j]>=abs(valAll[j])
		pValAll[j]=(sum(hold)+1)/(numpnts(hold)+1)
	endfor
	
	concatenate /NP=1 "valAll;",dimVals
	concatenate /NP=1 "pValAll;",pValSet
	
	duplicate /o pValSet xVals,yVals,hold
	xVals=p
	yVals=q
	redimension /N=(numpnts(pValSet)) hold,xVals,yVals
	
	sort hold,hold,xVals,yVals 				
	hold*=numpnts(hold)-p
	hold=limit(hold[p],0,1)
	
	for(i=0;i<numpnts(hold);i+=1)
		pValSet[xVals[i]][yVals[i]]=hold[i]
	endfor
end

function doBootAll(wv1,wv2,allWv)
	wave wv1,wv2,allWv
	make /o/n=(numpnts(allWv)) order=limit(floor(abs(enoise(numpnts(allWv)))),0,numpnts(allWv)-1)
	wv1=allWv[order[p]]
	wv2=allWv[order[p+numpnts(wv1)]]
end


function getMag(wv1,wv2,wv3)
	wave wv1,wv2,wv3
	variable val1=median(wv1)
	variable val2=median(wv2)
	variable val3=median(wv3)
	variable val=(val1-val2)+(val1-val3)
	return val
end


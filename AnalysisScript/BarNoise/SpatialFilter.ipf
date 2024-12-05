

pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

CONSTANT LinesSize=0.04, DeltaRF=.01,LengthRF=0.5,Step=10
CONSTANT PreFrame=15,PostFrame=15,StimFrame=1500,FrameTime=0.0333296423430406
CONSTANT ResponseSmpl=0.0001
CONSTANT SmoothTm=10

function doMany(howMany)
	variable howMany
	string pthStr="where are the stimulus and response files"
	newpath /M=pthStr/o/q filePath
	wave resp0
	if(waveExists(resp0))
		killwaves /z resp0
	endif
	loadWave /q/g/m/d/p=filePath/a=resp "resp.txt"
	wave resp0 
	setscale /p y,0,ResponseSmpl,resp0
	make /o/n=(dimsize(resp0,1)) linesRF_v
	setscale /p x,dimoffset(resp0,1),dimdelta(resp0,1),linesRF_v
	variable i
	for(i=0;i<howMany;i+=1)
		loadWave /q/g/m/d/p=filePath/n=stim "stim"+num2str(i)+".txt"
		wave stim0
		make /o/n=(dimsize(stim0,0),1,dimsize(stim0,1)) lines
		setscale /p z,0,FrameTime,lines
		lines[][0][]=stim0[p][r]
		killwaves /z stim0
		make /o/n=(dimsize(resp0,1)) linesRF_v
		setscale /p x,dimoffset(resp0,1),dimdelta(resp0,1),linesRF_v
		linesRF_v=resp0[i][p]
		if(deltax(linesRF_v)<0.001)
			reduceSamplingInPlace("linesRF_v",0.001)
		endif
		variable preTime=PreFrame*FrameTime
		deletepoints 0,preTime/deltax(linesRF_v),linesRF_v
		variable stimTime=StimFrame*FrameTime
		variable postTime=postFrame*FrameTime
		deletepoints stimTime/deltax(linesRF_v),numpnts(linesRF_v),linesRF_v
		getLinesRF()
		wave linesRF
		if(i==0)
			duplicate /o linesRF RF
			duplicate /o linesRF_v allResp
		else
			RF+=linesRF[p][q]
			concatenate /NP=1 "linesRF_v;",allResp
		endif
	endfor
	RF/=howMany
	
	make /o/n=(dimsize(RF,1),1,dimsize(RF,0)) holdRF
	setscale /p z,dimoffset(RF,0),dimdelta(RF,0),holdRF
	setscale /p x,dimoffset(RF,1),dimdelta(RF,1),holdRF
	holdRF=RF(z)[p]
	
	variable threshold=0
	variable totTime
	variable tmPnts
	for(i=0;i<howMany;i+=1)
		loadWave /q/g/m/d/p=filePath/n=stim "stim"+num2str(i)+".txt"
		wave stim0
		make /o/n=(dimsize(stim0,0),1,dimsize(stim0,1)) lines
		setscale /p z,0,FrameTime,lines
		lines[][0][]=stim0[p][r]
		killwaves /z stim0
		make /o/n=(dimsize(resp0,1)) linesRF_v
		setscale /p x,dimoffset(resp0,1),dimdelta(resp0,1),linesRF_v
		linesRF_v=resp0[i][p]
		if(deltax(linesRF_v)<0.001)
			reduceSamplingInPlace("linesRF_v",0.001)
		endif
		deletepoints 0,PreTime/deltax(linesRF_v),linesRF_v
		deletepoints StimTime/deltax(linesRF_v),PostTime/deltax(linesRF_v),linesRF_v
		if(dimsize(RF,1)>1)
			threshRF(holdRF,threshold,.025,.3)
			wave rfx,rfy
		else
			make /o/n=1 rfx=0,rfy=0
		endif
		totTime=rightx(linesRF_v)-deltax(linesRF_v)
		make /o/n=(dimsize(lines,0),dimsize(lines,1),floor(totTime/dimdelta(holdRF,2))) stm
		setscale  /p z,0,dimdelta(holdRF,2),stm
		stm=lines[p][q](z)
		wavestats /q stm
		stm-=v_avg
		stm/=v_avg
		convolverf (stm,holdRF,rfx,rfy)
		wave convsum
		tmPnts=(LengthRF/deltax(convsum))
		deletepoints 0,tmPnts, convsum
		setscale /p x,LengthRF,deltax(convsum),convsum
		duplicate /o convsum oneResp
		oneResp=linesRF_v(x)
		wavestats /q oneResp
		duplicate /o oneResp smHold
		smooth /b (SmoothTm/deltax(oneResp)),smhold
		oneResp-=smhold
		if(i==0)
			duplicate /o convsum lp
			duplicate /o oneResp v
			duplicate /o stm allStim
		else
			concatenate /NP=0 "convsum;",lp
			concatenate /NP=0 "oneResp;",v
			concatenate /NP=2 "stm;",allStim
		endif
	endfor
	wavestats /q allStim
	variable cont=v_sdev
	wavestats /q lp
	RF*=cont
	RF/=v_sdev
	lp*=cont
	lp/=v_sdev
	
	doNL(lp,v,50)
	wave NLy,NLx,NLySD
	duplicate /o NLy vNLy
	duplicate /o NLx vNLx
	duplicate /o NLySD vNLySD
	duplicate /o lp vlp
	
	RF/=dimdelta(RF,0)
	duplicate /o RF vRF
	
	killWindow /z linesGraphV
	doLinesGraphV()
	
	killwaves /z resp0,holdRF,stm,convsum,lp,v,rfx,rfy,rfmags,rf_thr,rf1p,smHld,st1p,vpmean,vpn,vpsd,vpsem,vpx,vpy,bend,bst,cfit,conv1p
end

function getlinesRF()
	wave linesRF_v
	wave lines
	wavestats /q lines 
	variable avgStim=v_avg
	make /o/n=(LengthRF/DeltaRF,dimsize(lines,0)) linesRF=0
	setscale /p x,0,DeltaRF,linesRF
	setscale /p y,0,LinesSize,linesRF
	duplicate /o linesRF_v hold,smHold
	smooth /b (SmoothTm/deltax(linesRF_v)),smhold
	hold-=smhold
	variable convert
	variable start=LengthRF/deltax(linesRF_v)
	variable i,j
	for(i=start;i<numpnts(linesRF_v);i+=step)
		convert=i*deltax(linesRF_v)
		linesRF+=(lines[q][0](convert-x)-avgStim)*hold[i]/avgStim
	endfor
end

Function ReduceSamplingInPlace(Nwavein, binSize)
	string Nwavein
	variable binSize

	wave wavein=$Nwavein
	Redimension/S wavein
	variable originalbin, reduction,  startbin
	originalbin=deltax(wavein)

	if (originalbin>=binSize)
		print "input wave is undersampled"
		return 0
	endif

	string Nwaveout
	sprintf Nwaveout, "%s_filtered", Nwavein
	make/o/n=(floor((numpnts(wavein)*originalbin)/binSize)) $Nwaveout
	wave waveout=$Nwaveout
	startbin=leftx($nwavein)+binSize/2

	SetScale/P x, startbin,binSize, waveout
	waveout=mean(wavein, x-binSize/2, x+binSize/2)
	duplicate /o waveout,$Nwavein
	killwaves waveout
end

function convolverf (wst,wrf,wrfx,wrfy)
	wave wst
	wave wrf
	wave wrfx
	wave wrfy
	variable pix
	//st1p:stimulus one pixel in milliseconds, number of points is
	//(number of frames in stimulus) * (frame time in seconds) * (1000 ms / s)
	//cor1p: convolution for one pixel, points are in ms
	//corsum: sum of all pixel 
	make /o/n=(dimsize(wst,2)) st1p,conv1p,convsum
	setscale /p x,0,dimdelta(wst,2),st1p,conv1p,convsum
	//rf1p:time course of one pixel of receptive field in ms, number of points is
	//(number of timepoints in receptive field) * (delta-t in seconds) * (1000 ms / s)
	make /o/n=(dimsize(wrf,2)) rf1p
	setscale /p x,0,dimdelta(wrf,2),rf1p
	convsum=0
	for (pix=0;pix<numpnts (wrfx);pix+=1) //Loop over number of significant pixels in RF
		rf1p=wrf[wrfx[pix]][wrfy[pix]](x)
		st1p=wst[wrfx[pix]][wrfy[pix]](x)
		duplicate /o st1p,conv1p
		convolve rf1p,conv1p
		convsum+=conv1p
		deletepoints 0,numpnts(rf1p),conv1p
	endfor
end

function threshrf (rfin,thr,tstart,tend)
	wave rfin
	variable thr
	variable tstart
	variable tend
	//rfmags:root-mean-squared value  of each pixel
	//rfx:wave containing x coord of significant pixels
	//rfy:wave containing y coord of significant pixels
	make /o/n=(dimsize (rfin,0)*dimsize(rfin,1)) rfmags,rfx,rfy
	variable px,py,pix
	make /o/n=(dimsize(rfin,2)) onepixtimesq // For one pixel p, holds p(t)^2, used to
	                                                                         // compute the rms value of that pixel over time
	setscale /p x,dimoffset(rfin,2),dimdelta(rfin,2),onepixtimesq
	for (px=0;px<dimsize(rfin,0);px+=1)        //Loop over x coord
		for (py=0;py<dimsize(rfin,1);py+=1) //Loop over y coord
			onepixtimesq=rfin[px][py][p]^2
			wavestats /q/r=(tstart,tend) onepixtimesq 
			rfmags[pix]=sqrt(v_avg) //rms value of one pix
			rfx[pix]=px
			rfy[pix]=py
			pix+=1
		endfor
	endfor
	sort /r rfmags,rfmags,rfx,rfy //sort the pixels with the greatest rms values
	wavestats /q  rfmags
	rfx=selectnumber( rfmags[p]>v_sdev*thr,nan,rfx) //get pixels greater than the thr in std devs
	rfy=selectnumber( rfmags[p]>v_sdev*thr,nan,rfy) //get pixels greater than the thr in std devs
	rfmags=selectnumber( rfmags[p]>v_sdev*thr,nan,rfmags)
	wavestats /q rfmags
	deletepoints v_npnts,numpnts(rfmags),rfmags,rfx,rfy 
	make /o/n=(dimsize(rfin,0),dimsize(rfin,1)) rf_thr //thresholded spatial receptive field
	rf_thr=0
	for (pix=0;pix<v_npnts;pix+=1)
		rf_thr[rfx[pix]][rfy[pix]]=rfmags[pix]
	endfor
end

function doNL(xwv,ywv,pnts)
	wave xwv,ywv
	variable pnts
	duplicate /o xwv xSort
	duplicate /o ywv ySort
	sort xSort,xSort,ySort
	make /o/n=(ceil(numpnts(xSort)/pnts)) NLx,NLy,NLySd
	make /o/n=(pnts) NLx,NLy,NLySd
	variable stepSz=ceil(numpnts(xSort)/pnts)
	variable i
	for(i=0;i<pnts;i+=1)
		wavestats /q/r=[i*stepSz,(i+1)*stepSz-1] xSort
		NLx[i]=v_avg
		wavestats /q/r=[i*stepSz,(i+1)*stepSz-1] ySort
		NLy[i]=v_avg
		NLySd[i]=v_sem
	endfor
	killwaves /z xSort,ySort
end

function doLN()
	wave allResp
	
	string pthStr="where is the file that you want to load?"
	newpath /M=pthStr/o/q thisPath
	pathInfo thisPath
	string filePath=s_path

	variable i,j
	for(i=0;i<dimsize(allResp,1);i+=1)
		wave thisSpike=$"st"+num2str(i)
		
		loadWave /q/g/m/d/n=stim filePath+"stim"+num2str(i)+".txt"
		wave stim0
		make /o/n=(dimsize(stim0,0),1,dimsize(stim0,1)) lines
		setscale /p z,0,FrameTime,lines
		lines[][0][]=stim0[p][r]
		killwaves /z stim0
		wavestats /q lines
		lines-=v_avg
		lines/=v_avg
		
		make /o/n=(LengthRF/DeltaRF,dimsize(lines,0)) spikeFilt=0
		setscale /p x,0,DeltaRF,spikeFilt
		setscale /p y,0,LinesSize,spikeFilt
		for(j=0;j<numpnts(thisSpike);j+=1)
			if(thisSpike[j]>LengthRF)
				spikeFilt+=lines[q][0](thisSpike[j]-x)
			endif
		endfor
		if(i==0)
			duplicate /o spikeFilt spKer
		else
			spKer+=spikeFilt
		endif
	endfor
	killwaves /z spikeFilt
	
	make /o/n=(dimsize(spKer,1),1,dimsize(spKer,0)) holdRF
	setscale /p x,dimoffset(spKer,1),dimdelta(spKer,1),holdRF
	setscale /p z,dimoffset(spKer,0),dimdelta(spKer,0),holdRF
	holdRF[][0][]=spKer[r][p]
	
	variable threshold=0
	variable totTime
	variable tmPnts
	for(i=0;i<dimsize(allResp,1);i+=1)
		loadWave /q/g/m/d/n=stim filePath+"stim"+num2str(i)+".txt"
		wave stim0
		make /o/n=(dimsize(stim0,0),1,dimsize(stim0,1)) lines
		setscale /p z,0,FrameTime,lines
		lines[][0][]=stim0[p][r]
		killwaves /z stim0
		if(dimsize(spKer,1)>1)
			threshRF(holdRF,threshold,.025,.3)
			wave rfx,rfy
		else
			make /o/n=1 rfx=0,rfy=0
		endif
		totTime=StimFrame*FrameTime
		make /o/n=(dimsize(lines,0),dimsize(lines,1),floor(totTime/dimdelta(holdRF,2))) stm
		setscale  /p z,0,dimdelta(holdRF,2),stm
		stm=lines[p][q](z)
		wavestats /q stm
		stm-=v_avg
		stm/=v_avg
		convolverf (stm,holdRF,rfx,rfy)
		wave convsum
		tmPnts=(LengthRF/deltax(convsum))
		deletepoints 0,tmPnts, convsum
		setscale /p x,LengthRF,deltax(convsum),convsum
		duplicate /o convsum oneResp
		wave thisSpike=$"st"+num2str(i)
		histogram /b=2 thisSpike,oneResp
		oneResp/=deltax(oneResp)
		if(i==0)
			duplicate /o convsum lp
			duplicate /o oneResp v
			duplicate /o stm allStim
		else
			concatenate /NP=0 "convsum;",lp
			concatenate /NP=0 "oneResp;",v
			concatenate /NP=2 "stm;",allStim
		endif
	endfor
	wavestats /q allStim
	variable cont=v_sdev
	wavestats /q lp
	spKer*=cont
	spKer/=v_sdev
	lp*=cont
	lp/=v_sdev
	
	doNL(lp,v,50)
	wave NLy,NLx,NLysd
	duplicate /o NLy spNLy
	duplicate /o NLySd spNLySd
	duplicate /o NLx spNLx
	duplicate /o lp spLP
	
	spKer/=dimdelta(spKer,0)
	duplicate /o spKer spRF
	
	killwaves /z holdRF,stm,convsum,lp,v,rfx,rfy,rfmags,rf_thr,rf1p,smHld,st1p,vpmean,vpn,vpsd,vpsem,vpx,vpy,bend,bst,cfit,conv1p
end

function extractSpikes()
	variable /g thr=15,d1thr=3,d2thr=3,fileNum=0
	wave allResp
	variable i
	for(i=0;i<dimsize(allResp,1);i+=1)
		make /o/n=(dimsize(allResp,0)) oneResp
		setscale /p x,dimoffset(allResp,0),dimdelta(allResp,0),oneResp
		oneResp=allResp[p][i]
		duplicate /o oneResp smResp
		smooth /b SmoothTm/deltax(oneResp),smResp
		oneResp-=smResp
		duplicate /o oneResp $"resp"+num2str(i)
	endfor
	killwindow /z spikeRemovalPanel
	doSpikeRemovalPanel()
	wave resp0
	spikeRemovalSetUp(resp0)
	killwaves /z oneResp,smResp
end

function doSpikeRemovalPanel()
	NVAR thr,d1thr,d2thr
	NVAR fileNum
	wave allResp
	NewPanel /k=1/n=spikeRemovalPanel/W=(350,45,596,211)
	Button b2,pos={130,102},size={104,20},proc=RemoveProc,title="Remove Spikes"
	Button b3,pos={155,135},size={54,20},proc=FinishProc,title="Finish"
	SetVariable sv0,pos={65,10},size={100,15},limits={0,dimsize(allResp,1)-1,1},fsize=14,proc=NextFile,value=fileNum,title="which file"
	SetVariable sv1,pos={10,43},size={145,15},limits={-80,80,0.5},fsize=14,proc=UpdateProc,value= thr,title="mem thresh"
	SetVariable sv2,pos={10,73},size={110,15},limits={0,10,0.5},fsize=14,proc=UpdateProc,value= d1thr,title="d1 thresh"
	SetVariable sv3,pos={10,103},size={110,15},limits={0,10,0.5},fsize=14,proc=UpdateProc,value= d2thr,title="d2 thresh"
End

Function UpdateProc(SV_Struct) : SetVariableControl
	STRUCT WMSetVariableAction &SV_Struct
	wave memThresh,d1Thresh,d2Thresh
	NVAR thr,d1thr,d2thr
	memThresh[1][]=thr
	d1thresh[2][]=-d1thr
	d1thresh[3][]=d1thr
	d2thresh[1][]=-d2thr
	d2thresh[2][]=d2thr
	return 0
End

Function NextFile(SV_Struct) : SetVariableControl
	STRUCT WMSetVariableAction &SV_Struct
	NVAR fileNum
	wave thisResp=$"resp"+num2str(fileNum)
	spikeRemovalSetUp(thisResp)
	return 0
End

Function RemoveProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	NVAR thr,d1thr,d2thr,fileNum
	wave thisResp=$"resp"+num2str(fileNum)
	switch( ba.eventCode )
		case 2: // mouse up
			removespikesint(thisResp,thr,d1thr,d2thr)
			wave wout,spiketimes, spNLy, spNLx, spRF
			duplicate /o wout $"wout"+num2str(fileNum)
			duplicate /o spiketimes $"st"+num2str(fileNum)
			duplicate /o spNLy $"spNLy"+num2str(fileNum)  // JL: 20210504 added to get NLy from each epoches
			duplicate /o spNLx $"spNLx"+num2str(fileNum)  // JL: 20210504 added to get NLy from each epoches
			duplicate /o spRF $"spRF"+num2str(fileNum)    // JL: 20210504 added to get NLy from each epoches
			break
		case -1: // control being killed
			break
	endswitch
	
	return 0
End

Function FinishProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	variable i
	wave allResp
	switch( ba.eventCode )
		case 2: // mouse up
			doLN()
			killWindow /z lineGraphSp
			doLinesGraphSp()
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

function spikeRemovalSetUp(win)
	wave win
	NVAR thr,d1thr,d2thr
	duplicate /o win,wraw,wdif,wdif2,wout
	differentiate wdif;wdif*=deltax(wdif)
	differentiate wdif2;differentiate wdif2;wdif2*=deltax(wdif)^2
	wavestats /q win
	wout=v_avg
	make /o/n=(4,2) d1thresh
	make /o/n=(3,2) d2thresh
	make /o/n=(2,2) memthresh
	wavestats /q wdif
	memthresh[0][0]=v_min
	memthresh[0][1]=v_max
	d2thresh[0][0]=v_min
	d2thresh[0][1]=v_max
	memthresh[1][]=thr
	wavestats /q wraw
	d1thresh[0][0]=v_min
	d1thresh[0][1]=v_max
	wavestats /q wdif2
	d1thresh[1][0]=v_min
	d1thresh[1][1]=v_max
	d1thresh[2][]=-d1thr
	d1thresh[3][]=d1thr
	d2thresh[1][]=-d2thr
	d2thresh[2][]=d2thr
	
	doWindow spikeRemovalGraph
	if(v_flag==0)
		doSpikeRemovalGraph()
	endif
end

//Removes spikes
//thr0:threshold for waves value
//thr1:threshold for absolute value of first derivative
//thr2:threhsold for absolute value of second derivative
//Points that are below all of these thresholds are collected/
//Then the wave is interpolated
function removespikesint(win,thr0,thr1,thr2)
	wave win
	variable thr0,thr1,thr2
	duplicate /o win,wselect,wdif,wdif2,wspike
	differentiate wdif;wdif*=deltax(wdif)
	differentiate wdif2;differentiate wdif2;wdif2*=deltax(wdif2)^2
	wselect=selectnumber((win<thr0)%&(abs(wdif)<thr1)%&(abs(wdif2)<thr2),nan,win)
	wspike=selectnumber((win<thr0)%&(abs(wdif)<thr1)%&(abs(wdif2)<thr2),1,0)
	findlevels /q/EDGE=1 wspike,1
	wave w_findlevels
	duplicate /o w_findlevels spikeTimes
	interpolate2 /t=3/f=0.01/y=wout wselect
end

Function doSpikeRemovalGraph()
	wave wraw,wout,wdif,wdif2,memThresh,d1Thresh,d2Thresh
	Display /k=1/n=spikeRemovalGraph/W=(600,45,1126,229) wdif vs wraw
	AppendToGraph/L=l2/B=b2 wraw,wout
	AppendToGraph/B=b1 wdif vs wdif2
	AppendToGraph memthresh[0][*] vs memthresh[1][*]
	AppendToGraph d1thresh[2][*] vs d1thresh[0][*]
	AppendToGraph d1thresh[3][*] vs d1thresh[0][*]
	AppendToGraph/B=b1 d1thresh[2][*] vs d1thresh[1][*]
	AppendToGraph/B=b1 d1thresh[3][*] vs d1thresh[1][*]
	AppendToGraph/B=b1 d2thresh[0][*] vs d2thresh[1][*]
	AppendToGraph/B=b1 d2thresh[0][*] vs d2thresh[2][*]
	ModifyGraph mode(wdif)=2,mode(wdif#1)=2
	ModifyGraph rgb(wraw)=(0,0,0),rgb(wout)=(65535,43690,0),rgb(memthresh)=(0,0,0),rgb(d1thresh)=(0,0,0)
	ModifyGraph rgb(d1thresh#1)=(0,0,0),rgb(d1thresh#2)=(0,0,0),rgb(d1thresh#3)=(0,0,0)
	ModifyGraph rgb(d2thresh)=(0,0,0),rgb(d2thresh#1)=(0,0,0),rgb(wdif)=(1,52428,52428)
	ModifyGraph rgb(wdif#1)=(1,52428,52428)
	ModifyGraph nticks=3
	ModifyGraph fSize=10
	ModifyGraph standoff=0
	ModifyGraph lblPosMode(bottom)=1,lblPosMode(b1)=1,lblPosMode(l2)=3,lblPosMode(b2)=1
	ModifyGraph lblPos(left)=32,lblPos(bottom)=27,lblPos(l2)=35
	ModifyGraph ZisZ=1
	ModifyGraph zapTZ=1
	ModifyGraph btLen=1.5
	ModifyGraph freePos(b1)={0,kwfraction}
	ModifyGraph freePos(l2)={0.67,kwFraction}
	ModifyGraph freePos(b2)={0,kwfraction}
	ModifyGraph axisEnab(bottom)={0,0.25}
	ModifyGraph axisEnab(b1)={0.3,0.55}
	ModifyGraph axisEnab(b2)={0.67,1}
	Label left "\\Z12\\F'Helvetica'dmV/dt"
	Label bottom "\\Z12\\F'Helvetica'mV"
	Label b1 "\\Z12\\F'Helvetica'd\\S2\\M\\Z12\\F'Helvetica'mV/dt\\S2"
	Label l2 "\\Z12\\F'Helvetica'mV"
	Label b2 "\\Z12\\F'Helvetica'Time (s)"
End

function doLinesGraphV()
	wave vRF,vNLy,vNLx,vNLysd
	Display /k=0/n=linesGraphV/W=(12,175,259,270)/L=l1/B=b1 vNLy vs vNLx
	AppendImage vRF
	ModifyImage vRF ctab= {-5,5,RedWhiteBlue,0}
	ModifyGraph mode=3
	ModifyGraph marker=19
	ModifyGraph rgb=(0,0,0)
	ModifyGraph mirror(left)=0,mirror(bottom)=0
	ModifyGraph nticks(left)=2,nticks(bottom)=2,nticks(l1)=2,nticks(b1)=4
	ModifyGraph fSize=10
	ModifyGraph standoff=0
	ModifyGraph lblPosMode(bottom)=1,lblPosMode(l1)=3,lblPosMode(b1)=1
	ModifyGraph lblPos(left)=42,lblPos(bottom)=25,lblPos(l1)=28
	ModifyGraph ZisZ=1
	ModifyGraph zapTZ=1
	ModifyGraph btLen=1.5
	ModifyGraph freePos(l1)={0.65,kwFraction}
	ModifyGraph freePos(b1)={0,kwFraction}
	ModifyGraph axisEnab(bottom)={0,0.45}
	ModifyGraph axisEnab(b1)={0.65,1}
	Label left "\\Z12\\F'Helvetica'Size (mm)"
	Label bottom "\\Z12\\F'Helvetica'Delay (s)"
	Label l1 "\\Z12\\F'Helvetica'Output (mV)"
	Label b1 "\\Z12\\F'Helvetica'Input"
	ErrorBars/T=0/L=2 vNLy Y,wave=(vNLySd,vNLySd)
End

function doLinesGraphSp()
	wave spRF,spNLy,spNLx,spNLysd
	Display /k=1/n=linesGraphSp/W=(14,296,261,391)/L=l1/B=b1 spNLy vs spNLx
	AppendImage spRF
	ModifyImage spRF ctab= {-5,5,RedWhiteBlue,0}
	ModifyGraph mode=3
	ModifyGraph marker=19
	ModifyGraph rgb=(0,0,0)
	ModifyGraph mirror(left)=0,mirror(bottom)=0
	ModifyGraph nticks(left)=2,nticks(bottom)=2,nticks(l1)=2,nticks(b1)=4
	ModifyGraph fSize=10
	ModifyGraph standoff=0
	ModifyGraph lblPosMode(bottom)=1,lblPosMode(l1)=3,lblPosMode(b1)=1
	ModifyGraph lblPos(left)=42,lblPos(bottom)=25,lblPos(l1)=28
	ModifyGraph ZisZ=1
	ModifyGraph zapTZ=1
	ModifyGraph btLen=1.5
	ModifyGraph freePos(l1)={0.65,kwFraction}
	ModifyGraph freePos(b1)={0,kwFraction}
	ModifyGraph axisEnab(bottom)={0,0.45}
	ModifyGraph axisEnab(b1)={0.65,1}
	Label left "\\Z12\\F'Helvetica'Size (mm)"
	Label bottom "\\Z12\\F'Helvetica'Delay (s)"
	Label l1 "\\Z12\\F'Helvetica'Output (mV)"
	Label b1 "\\Z12\\F'Helvetica'Input"
	ErrorBars/T=0/L=2 spNLy Y,wave=(spNLySd,spNLySd)
End

function doGraphSp()
	wave spRF,spNLy,spNLx,spNLysd
	Display /k=1/n=linesGraphSp/W=(14,296,261,391)/L=l1/B=b1 spNLy vs spNLx
	AppendImage spRF
	ModifyImage spRF ctab= {-5,5,RedWhiteBlue,0}
	ModifyGraph mode=3
	ModifyGraph marker=19
	ModifyGraph rgb=(0,0,0)
	ModifyGraph mirror(left)=0,mirror(bottom)=0
	ModifyGraph nticks(left)=2,nticks(bottom)=2,nticks(l1)=2,nticks(b1)=4
	ModifyGraph fSize=10
	ModifyGraph standoff=0
	ModifyGraph lblPosMode(bottom)=1,lblPosMode(l1)=3,lblPosMode(b1)=1
	ModifyGraph lblPos(left)=42,lblPos(bottom)=25,lblPos(l1)=28
	ModifyGraph ZisZ=1
	ModifyGraph zapTZ=1
	ModifyGraph btLen=1.5
	ModifyGraph freePos(l1)={0.65,kwFraction}
	ModifyGraph freePos(b1)={0,kwFraction}
	ModifyGraph axisEnab(bottom)={0,0.45}
	ModifyGraph axisEnab(b1)={0.65,1}
	Label left "\\Z12\\F'Helvetica'Size (mm)"
	Label bottom "\\Z12\\F'Helvetica'Delay (s)"
	Label l1 "\\Z12\\F'Helvetica'Output (mV)"
	Label b1 "\\Z12\\F'Helvetica'Input"
	ErrorBars/T=0/L=2 spNLy Y,wave=(spNLySd,spNLySd)
End
#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Written by David B. Kastner
// Editd by Joo Yeun Lee

function getSpaceSP(wv)  //getspace spikes
	wave wv
	wv*=deltax(wv)
	matrixMultiply wv,wv/t
	wave m_product
	matrixEigenV /SYM/EVEC m_product
	wave m_eigenvectors,w_eigenvalues
	duplicate /o w_eigenvalues var
	wavestats /q var
	var/=v_sum
	make /o/n=(dimsize(m_eigenvectors,0)) spPC1,spPC2,spPC3
	setscale /p x,dimoffset(wv,0),dimdelta(wv,0),spPC1,spPC2,spPC3
	spPC1=m_eigenvectors[p][dimsize(m_eigenvectors,1)-1]
	spPC2=m_eigenvectors[p][dimsize(m_eigenvectors,1)-2]
	spPC3=m_eigenvectors[p][dimsize(m_eigenvectors,1)-3]
	matrixMultiply wv/t,spPC1
	make /o/n=(dimsize(m_product,0)) spspace1,spspace2,spspace3
	setscale /p x,dimoffset(wv,1),dimdelta(wv,1),spspace1,spspace2,spspace3
	spspace1=m_product[p][0]
	matrixMultiply wv/t,spPC2
	spspace2=m_product[p][0]
	matrixMultiply wv/t,spPC3
	spspace3=m_product[p][0]
	wavestats /q spspace1
	wv/=deltax(wv)
	
	killWindow /z PCgraph
	doPCgraphSP()
end

function doPCgraphSP()
	wave spPC1,spPC2,spPC3
	wave spspace1,spspace2,spspace3
	duplicate /o spPC1 zeroPC
	zeroPC=0
	duplicate /o spspace1 zeroSpace
	zeroSpace=0
	wave var
	wavestats /q spPC1
	variable PCval=max(v_max,abs(v_min))
	wavestats /q spPC2
	PCval=max(v_max,PCval)
	PCval=max(PCval,abs(v_min))
	wavestats /q spPC3
	PCval=max(v_max,PCval)
	PCval=max(PCval,abs(v_min))
	wavestats /q spspace1
	variable spMax=max(0,v_max)
	variable spMin=min(0,v_min)
	wavestats /q spspace2
	spMax=max(spMax,v_max)
	spMin=min(spMin,v_min)
	wavestats /q spspace3
	spMax=max(spMax,v_max)
	spMin=min(spMin,v_min)
	Display /k=1/n=PCgraph/W=(364,512,701,768) zeroPC,spPC3
	AppendToGraph/L=l1 zeroPC,spPC2
	AppendToGraph/L=l2 zeroPC,spPC1
	AppendToGraph/L=l3/B=b1 zeroSpace,spspace3
	AppendToGraph/L=l4/B=b1 zeroSpace,spspace2
	AppendToGraph/L=l5/B=b1 zeroSpace,spspace1
	ModifyGraph rgb=(0,0,0)
	ModifyGraph nticks=2
	ModifyGraph lowTrip=0.01
	ModifyGraph lblMargin(left)=5,lblMargin(l1)=5,lblMargin(l2)=5
	ModifyGraph standoff=0
	ModifyGraph lblPosMode(left)=1,lblPosMode(bottom)=1,lblPosMode(l1)=1,lblPosMode(l2)=1
	ModifyGraph lblPosMode(b1)=1
	ModifyGraph lblPos(left)=60,lblPos(bottom)=27
	ModifyGraph lstyle(zeroPC)=1,lstyle(zeroPC#1)=1,lstyle(zeroPC#2)=1,lstyle(zeroSpace)=1,lstyle(zeroSpace#1)=1,lstyle(zeroSpace#2)=1
	ModifyGraph ZisZ=1
	ModifyGraph zapTZ=1
	ModifyGraph btLen=1.5
	ModifyGraph freePos(l1)={0,kwfraction}
	ModifyGraph freePos(l2)={0,kwfraction}
	ModifyGraph freePos(l3)={0.6,kwFraction}
	ModifyGraph freePos(b1)={0,kwfraction}
	ModifyGraph freePos(l4)={0.6,kwFraction}
	ModifyGraph freePos(l5)={0.6,kwFraction}
	ModifyGraph axisEnab(left)={0,0.3}
	ModifyGraph axisEnab(bottom)={0,0.4}
	ModifyGraph axisEnab(l1)={0.35,0.65}
	ModifyGraph axisEnab(l2)={0.7,1}
	ModifyGraph axisEnab(l3)={0,0.3}
	ModifyGraph axisEnab(b1)={0.6,1}
	ModifyGraph axisEnab(l4)={0.35,0.65}
	ModifyGraph axisEnab(l5)={0.7,1}
	Label left "\\Z12\\F'Helvetica'spPC3 "+num2str(round(var[numpnts(var)-3]*100))+"%"
	Label bottom "\\Z12\\F'Helvetica'Delay (s)"
	Label l1 "\\Z12\\F'Helvetica'spPC2 "+num2str(round(var[numpnts(var)-2]*100))+"%"
	Label l2 "\\Z12\\F'Helvetica'spPC1 "+num2str(round(var[numpnts(var)-1]*100))+"%"
	Label b1 "\\Z12\\F'Helvetica'Distance (mm)"
	SetAxis left -PCval,PCval
	SetAxis l1 -PCval,PCval
	SetAxis l2 -PCval,PCval
	SetAxis l3 spMin,spMax
	SetAxis l4 spMin,spMax
	SetAxis l5 spMin,spMax
End

function getSpaceMem(wv)   //getspace voltage
	wave wv
	wv*=deltax(wv)
	matrixMultiply wv,wv/t
	wave m_product
	matrixEigenV /SYM/EVEC m_product
	wave m_eigenvectors,w_eigenvalues
	duplicate /o w_eigenvalues var
	wavestats /q var
	var/=v_sum
	make /o/n=(dimsize(m_eigenvectors,0)) memPC1,memPC2,memPC3
	setscale /p x,dimoffset(wv,0),dimdelta(wv,0),memPC1,memPC2,memPC3
	memPC1=m_eigenvectors[p][dimsize(m_eigenvectors,1)-1]
	memPC2=m_eigenvectors[p][dimsize(m_eigenvectors,1)-2]
	memPC3=m_eigenvectors[p][dimsize(m_eigenvectors,1)-3]
	matrixMultiply wv/t,memPC1
	make /o/n=(dimsize(m_product,0)) memspace1,memspace2,memspace3
	setscale /p x,dimoffset(wv,1),dimdelta(wv,1),memspace1,memspace2,memspace3
	memspace1=m_product[p][0]
	matrixMultiply wv/t,memPC2
	memspace2=m_product[p][0]
	matrixMultiply wv/t,memPC3
	memspace3=m_product[p][0]
	wavestats /q memspace1
	wv/=deltax(wv)
	
	killWindow /z PCgraph
	doPCgraphMem()
end

function doPCgraphMem()
	wave memPC1,memPC2,memPC3
	wave memspace1,memspace2,memspace3
	duplicate /o memPC1 zeroPC
	zeroPC=0
	duplicate /o memspace1 zeroSpace
	zeroSpace=0
	wave var
	wavestats /q memPC1
	variable PCval=max(v_max,abs(v_min))
	wavestats /q memPC2
	PCval=max(v_max,PCval)
	PCval=max(PCval,abs(v_min))
	wavestats /q memPC3
	PCval=max(v_max,PCval)
	PCval=max(PCval,abs(v_min))
	wavestats /q memspace1
	variable spMax=max(0,v_max)
	variable spMin=min(0,v_min)
	wavestats /q memspace2
	spMax=max(spMax,v_max)
	spMin=min(spMin,v_min)
	wavestats /q memspace3
	spMax=max(spMax,v_max)
	spMin=min(spMin,v_min)
	Display /k=1/n=PCgraph/W=(364,512,701,768) zeroPC,memPC3
	AppendToGraph/L=l1 zeroPC,memPC2
	AppendToGraph/L=l2 zeroPC,memPC1
	AppendToGraph/L=l3/B=b1 zeroSpace,memspace3
	AppendToGraph/L=l4/B=b1 zeroSpace,memspace2
	AppendToGraph/L=l5/B=b1 zeroSpace,memspace1
	ModifyGraph rgb=(0,0,0)
	ModifyGraph nticks=2
	ModifyGraph lowTrip=0.01
	ModifyGraph lblMargin(left)=5,lblMargin(l1)=5,lblMargin(l2)=5
	ModifyGraph standoff=0
	ModifyGraph lblPosMode(left)=1,lblPosMode(bottom)=1,lblPosMode(l1)=1,lblPosMode(l2)=1
	ModifyGraph lblPosMode(b1)=1
	ModifyGraph lblPos(left)=60,lblPos(bottom)=27
	ModifyGraph lstyle(zeroPC)=1,lstyle(zeroPC#1)=1,lstyle(zeroPC#2)=1,lstyle(zeroSpace)=1,lstyle(zeroSpace#1)=1,lstyle(zeroSpace#2)=1
	ModifyGraph ZisZ=1
	ModifyGraph zapTZ=1
	ModifyGraph btLen=1.5
	ModifyGraph freePos(l1)={0,kwfraction}
	ModifyGraph freePos(l2)={0,kwfraction}
	ModifyGraph freePos(l3)={0.6,kwFraction}
	ModifyGraph freePos(b1)={0,kwfraction}
	ModifyGraph freePos(l4)={0.6,kwFraction}
	ModifyGraph freePos(l5)={0.6,kwFraction}
	ModifyGraph axisEnab(left)={0,0.3}
	ModifyGraph axisEnab(bottom)={0,0.4}
	ModifyGraph axisEnab(l1)={0.35,0.65}
	ModifyGraph axisEnab(l2)={0.7,1}
	ModifyGraph axisEnab(l3)={0,0.3}
	ModifyGraph axisEnab(b1)={0.6,1}
	ModifyGraph axisEnab(l4)={0.35,0.65}
	ModifyGraph axisEnab(l5)={0.7,1}
	Label left "\\Z12\\F'Helvetica'memPC3 "+num2str(round(var[numpnts(var)-3]*100))+"%"
	Label bottom "\\Z12\\F'Helvetica'Delay (s)"
	Label l1 "\\Z12\\F'Helvetica'memPC2 "+num2str(round(var[numpnts(var)-2]*100))+"%"
	Label l2 "\\Z12\\F'Helvetica'memPC1 "+num2str(round(var[numpnts(var)-1]*100))+"%"
	Label b1 "\\Z12\\F'Helvetica'Distance (mm)"
	SetAxis left -PCval,PCval
	SetAxis l1 -PCval,PCval
	SetAxis l2 -PCval,PCval
	SetAxis l3 spMin,spMax
	SetAxis l4 spMin,spMax
	SetAxis l5 spMin,spMax
End
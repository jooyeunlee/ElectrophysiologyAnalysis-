#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Written by David B. Kastner

function getSpace(wv)
	wave wv
	wv*=deltax(wv)
	matrixMultiply wv,wv/t
	wave m_product
	matrixEigenV /SYM/EVEC m_product
	wave m_eigenvectors,w_eigenvalues
	duplicate /o w_eigenvalues var
	wavestats /q var
	var/=v_sum
	make /o/n=(dimsize(m_eigenvectors,0)) PC1,PC2,PC3
	setscale /p x,dimoffset(wv,0),dimdelta(wv,0),PC1,PC2,PC3
	PC1=m_eigenvectors[p][dimsize(m_eigenvectors,1)-1]
	PC2=m_eigenvectors[p][dimsize(m_eigenvectors,1)-2]
	PC3=m_eigenvectors[p][dimsize(m_eigenvectors,1)-3]
	matrixMultiply wv/t,PC1
	make /o/n=(dimsize(m_product,0)) space1,space2,space3
	setscale /p x,dimoffset(wv,1),dimdelta(wv,1),space1,space2,space3
	space1=m_product[p][0]
	matrixMultiply wv/t,PC2
	space2=m_product[p][0]
	matrixMultiply wv/t,PC3
	space3=m_product[p][0]
	wavestats /q space1
	wv/=deltax(wv)
	
	killWindow /z PCgraph
	doPCgraph()
end

function doPCgraph()
	wave PC1,PC2,PC3
	wave space1,space2,space3
	duplicate /o PC1 zeroPC
	zeroPC=0
	duplicate /o space1 zeroSpace
	zeroSpace=0
	wave var
	wavestats /q PC1
	variable PCval=max(v_max,abs(v_min))
	wavestats /q PC2
	PCval=max(v_max,PCval)
	PCval=max(PCval,abs(v_min))
	wavestats /q PC3
	PCval=max(v_max,PCval)
	PCval=max(PCval,abs(v_min))
	wavestats /q space1
	variable spMax=max(0,v_max)
	variable spMin=min(0,v_min)
	wavestats /q space2
	spMax=max(spMax,v_max)
	spMin=min(spMin,v_min)
	wavestats /q space3
	spMax=max(spMax,v_max)
	spMin=min(spMin,v_min)
	Display /k=1/n=PCgraph/W=(364,512,701,768) zeroPC,PC3
	AppendToGraph/L=l1 zeroPC,PC2
	AppendToGraph/L=l2 zeroPC,PC1
	AppendToGraph/L=l3/B=b1 zeroSpace,space3
	AppendToGraph/L=l4/B=b1 zeroSpace,space2
	AppendToGraph/L=l5/B=b1 zeroSpace,space1
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
	Label left "\\Z12\\F'Helvetica'PC3 "+num2str(round(var[numpnts(var)-3]*100))+"%"
	Label bottom "\\Z12\\F'Helvetica'Delay (s)"
	Label l1 "\\Z12\\F'Helvetica'PC2 "+num2str(round(var[numpnts(var)-2]*100))+"%"
	Label l2 "\\Z12\\F'Helvetica'PC1 "+num2str(round(var[numpnts(var)-1]*100))+"%"
	Label b1 "\\Z12\\F'Helvetica'Distance (mm)"
	SetAxis left -PCval,PCval
	SetAxis l1 -PCval,PCval
	SetAxis l2 -PCval,PCval
	SetAxis l3 spMin,spMax
	SetAxis l4 spMin,spMax
	SetAxis l5 spMin,spMax
End
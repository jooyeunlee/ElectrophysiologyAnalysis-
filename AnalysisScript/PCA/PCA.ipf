#pragma TextEncoding = "MacRoman"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Scatter Dot Plot>

function getPCA(wv)
	wave wv 
	matrixtranspose wv								
	make /o /n=(dimsize(wv,1),dimsize(wv,0)) wv_s=wv[q][p]
	sumdimension /D=0 wv
	wave w_sumdimension
	w_sumdimension/=dimsize(wv,0) 				
	wv_s-=w_sumdimension[p]						
	sumdimension /D=1 wv_s
	matrixMultiply wv_s,wv_s/t 					
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
	matrixMultiply wv_s/t,PC1
	make /o/n=(dimsize(m_product,0)) proj1, proj2
	setscale /p x,dimoffset(wv_s,1),dimdelta(wv,1),proj1, proj2
	proj1=m_product[p][0]
	matrixMultiply wv_s/t,PC2
	proj2=m_product[p][0]
end
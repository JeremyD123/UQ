[n,w] = nwspgr('KPN',1,16);

sum(w.*hermiteN(15,n).*hermiteN(15,n))
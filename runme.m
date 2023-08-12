load fncdata
fncno = 2 ;  % 1: kowalik  ;  2 :rosenbrock
fnc = string(fnclist(fncno)) ; 
low = data(fncno,1) ; 
up = data(fncno,2) ; 
dim = data(fncno,3) ; 
% out  =  CSA(fnc , mydata , N , D , low , up , Epk)
out  =  CSA(fnc , [ ] , 20 , dim , low , up , 200e3)

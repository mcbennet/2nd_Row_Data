***,Calculation for Be atom, singlet and triplet
memory,512,m
gthresh,twoint=1.0E-12
geometry={
1
Na
Na  0.0 0.0 0.0
}

basis={
ecp,Na,10,2,0
3;
 1  ,  4.31167777 ,   1.00000000 
 3  ,  1.92568851 ,   4.31167777 
 2  ,  1.54949836 ,  -2.08313670 
2;
 2  ,  5.37766645 ,   6.23406351 
 2  ,  1.40841434 ,   9.07593069 
2;
 2  ,  1.37994865 ,   3.23272398 
 2  ,  0.86245288 ,   2.49407937 
include,aug-cc-pV5Z.basis
}

include,states.proc

do i=1,2
if (i.eq.1) then
    GS
else if (i.eq.2) then
    EA
endif
scf(i)=energy
_CC_NORM_MAX=2.0
{uccsd(t),maxit=400;core}
ccsd(i)=energy
enddo

table,scf,ccsd
save
type,csv

***,Calculation for Be atom, singlet and triplet
memory,512,m
gthresh,twoint=1.0E-12
geometry={
1
Mg
Mg  0.0 0.0 0.0
}

basis={
ecp,Mg,10,2,0
3
  1  ,  6.04853770 ,   2.00000000 
  3  ,  2.79698905 ,  12.09707540 
  2  ,  2.54740820 , -17.10831330 
2
  2  ,  5.93601730 ,   6.42863103 
  2  ,  1.59289065 ,  14.19549140 
2
  2  ,  1.58396886 ,   3.31506900 
  2  ,  1.07729691 ,   4.40302510 
include,aug-cc-pVQZ.basis
}

include,states.proc

do i=1,2
if (i.eq.1) then
    GS
else if (i.eq.2) then
    IP1
endif
scf(i)=energy
_CC_NORM_MAX=2.0
{uccsd(t),maxit=400;core}
ccsd(i)=energy
enddo

table,scf,ccsd
save
type,csv

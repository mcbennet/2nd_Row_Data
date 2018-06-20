***,Calculation for Be atom, singlet and triplet
memory,512,m
gthresh,twoint=1.0E-12
geometry={
1
Si
Si  0.0 0.0 0.0
}

basis={
ecp,Si,10,2,0
3
  1  ,  3.54632218 ,   4.00000000 
  3  ,  4.10257346 ,  14.18528871 
  2  ,  3.03189783 , -14.55031990 
2
  2  ,  7.40540967 ,   6.03326125 
  2  ,  2.38741220 ,  23.08201380 
2
  2  ,  2.17523625 ,   3.34088948 
  2  ,  1.94078734 ,  10.10971550 
include,aug-cc-pV5Z.basis
}

include,states.proc

do i=1,5
if (i.eq.1) then
    GS
else if (i.eq.2) then
    IP1
else if (i.eq.3) then
    IP2
else if (i.eq.4) then
    IP3
else if (i.eq.5) then
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

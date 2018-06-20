***,Calculation for Be atom, singlet and triplet
memory,512,m


! s  0 x^2,y^2,z^2
! p  1 x
! p -1 y
! d -2 xy
! p  0 z
! d  1 xz
! d -1 yz
! f -2 xyz

print,orbitals=2

gparam,NOBUFF
gthresh,twoint=1.e-12
geometry={
1
Cl
Cl  0.0 0.0 0.0
}

basis={
ecp,Cl,10,2,0;
3;
1,22.716551732563,7.0;
3,78.571856852001,159.015862127941;
2,7.473524361828001,-15.6531065;
2;
2,17.237085734082,6.50888648;
2,4.311484472548,46.763467;
2;
2,11.382757039936,2.9946477;
2,3.832187618897,28.0170341;

include,aug-cc-pCVTZ.basis
}

include,ppstates.proc

do i=1,8
    if (i.eq.1) then
        N1
    else if (i.eq.2) then
        N2
    else if (i.eq.3) then
        N3
    else if (i.eq.4) then
        N4
    else if (i.eq.5) then
        N5
    else if (i.eq.6) then
        N6
    else if (i.eq.7) then
        N7
    else if (i.eq.8) then
        N8
    endif
    scf(i)=energy
    _CC_NORM_MAX=2.0
    {uccsd(t),maxit=100;core}
    ccsd(i)=energy
enddo

table,scf,ccsd
save
type,csv

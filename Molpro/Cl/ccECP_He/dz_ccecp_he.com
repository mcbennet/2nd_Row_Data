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
5;
1,0.0,-15.0;
1,0.0,7.0;
1,20.01071281,15.0;
3,21.44804874,300.16069215;
2,21.17991361,-151.85882234;
2;
2,44.02974125,202.03781785;
2,65.80681305,25.00902536;
1;
2,1.0,0.0;


include,aug-cc-pCVDZ.basis
}

include,ppstates_sc.proc

do i=0,8
    if (i.eq.0) then
        Nn
    else if (i.eq.1) then
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
    scf(i+1)=energy
    _CC_NORM_MAX=2.0
    {uccsd(t),maxit=100;core}
    ccsd(i+1)=energy
enddo

table,scf,ccsd
save
type,csv

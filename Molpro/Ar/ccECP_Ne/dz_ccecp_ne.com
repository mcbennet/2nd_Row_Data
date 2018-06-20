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
Ar
Ar  0.0 0.0 0.0
}

basis={
ecp,Ar,10,2,0;
3;
1,8.98034652,8.0;
3,26.23635698,71.84277216;
2,7.10917092,-22.1300043;
2;
2,26.06800734,7.67833187;
2,5.17516183,60.3130371;
2;
2,17.09350095,3.47705723;
2,4.64108131,36.4723622;


include,aug-cc-pCVDZ.basis
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

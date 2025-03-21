

HeegHy(E,D)={

N = E[16][6][1];
pN = factor(N)[,1];
res=0;
for(i=1,#pN,res = res+kronecker(D,pN[i]^factor(N)[i,2]));

return(res-#pN);
}

P0(E,D,P)={
HeegnerPoints(E,D);
z = ellpointtoz(E,P);
print(z/zsums[1]);
}


Listquadform(D)={

Cl = quadclassunit(D)[1];
L = Vec(0,Cl);
M = Vec(0,Cl);
P = Vec(0,Cl);

for(i=1, Cl, L[i]=qfbpowraw(quadclassunit(D)[3][1],i));
for(i=1, Cl, M[i]=Mat(L[i]));
for(i=1, Cl, P[i]=subst(Pol(L[i]),x,x/y)*y^2);

return([L,M,P]);

}

Listquadform1(D)={

Cl = quadclassunit(D)[1];
L = Vec(0,Cl);
M = Vec(0,Cl);
P = Vec(0,Cl);

L[1]=Vec(quadpoly(D));
M[1]=[L[1][1], L[1][2]/2;L[1][2]/2, L[1][3]];
P[1]=subst(quadpoly(D),x,x/y)*y^2;

return([L,M,P]);

}



sRep(P,j)={

S = Vec(0,#pN);
for(i=1,#pN,S[i]=Polsol(P,pN[i]^factor(N)[i,2])[j[i]+1]);
return(S);

}


Polsol(Q,q)={

s=List([]);

for(k=0,q-1, for(l=0,q-1,if(!Mod(subst(subst(Q,y,k),x,l),q) && (k+l),listput(s,[l,k]))));

for(n=2, q-1, for(ind=1, #s, if(ind<= #s,if(Mod(s[1]*n,q)==Mod(s[ind],q), listpop(s,ind) ))));
for(n=2, q-1, for(ind=1, #s, if(ind<= #s,if(Mod(s[2]*n,q)==Mod(s[ind],q), listpop(s,ind) ))));

return(s);

}

Polsol1(Q,q)={

s=List([]);

for(k=0,q-1, for(l=0,q-1,if(!Mod(subst(subst(Q,y,k),x,l),q) && (k+l),listput(s,Mod([l,k],q)))));

for(n=2, q-1, for(ind=1, #s, if(ind<= #s,if(Mod(s[1]*n,q)==Mod(s[ind],q), listpop(s,ind) ))));
for(n=2, q-1, for(ind=1, #s, if(ind<= #s,if(Mod(s[2]*n,q)==Mod(s[ind],q), listpop(s,ind) ))));

return(s);

}

Mattos(L)={

v11=lift(chinese(Mod(L[1][1]*L[2][1],pN[1]^factor(N)[1,2]),Mod(L[2][1]*L[1][1],pN[2]^factor(N)[2,2])));
v21=lift(chinese(Mod(L[1][2]*L[2][1],pN[1]^factor(N)[1,2]),Mod(L[2][2]*L[1][1],pN[2]^factor(N)[2,2])));

v12=-bezout(v11,v21)[2];
v22=bezout(v11,v21)[1];

M = [v11,v12;v21,v22];

return(M);

}


binaryvector(len,num)={

vtr = Vec(0,len);

for(i=1,#binary(num),vtr[i]=binary(num)[#binary(num)-i+1]);

return(vtr);

}

subf1(q)={                              
tot = 0;                                
for(i=1,#V-1,tot+= q^i*V[i]);         
return(tot);
}

Vdef(E1)={

V = Vec(0,800000);
for(i=1,#V-1,V[i]=ellak(E1,i)/i);

}

HeegnerPoints(E,D)={

N = E[16][6][1];
pN = factor(N)[,1];
qfD=0;

if(quadclassunit(D)[1]!=1,qfD = Listquadform(D), qfD = Listquadform1(D));

polForm = qfD[3];

svec = Vec(0,#polForm);
tauvec = Vec(0, #polForm);
qtauvec = Vec(0,#polForm);
zvec = Vec(0,#polForm);
Pvec = Vec(0,#polForm);
zsums = Vec(0,#polForm);
Psums = Vec(0,#polForm);
zsrep = Vec(0,2^#pN);

for(i=1,#polForm, svec[i]=listrep(E,D,i));
for(i=1,#polForm, tauvec[i]=svec[i][2]);
for(i=1,#polForm, qtauvec[i]=svec[i][3]);
for(i=1,#polForm, zvec[i]=svec[i][4]);
for(i=1,#polForm, Pvec[i]=svec[i][5]);
for(i=1,#polForm, zsums[i]=svec[i][6]);
for(i=1,#polForm, Psums[i]=svec[i][7]);

for(i=1, 2^#pN, for(j=1,#polForm, zsrep[i]+=zvec[j][i]));

return(svec);

}


HeegnerPoints2(E,D)={

N = E[16][6][1];
pN = factor(N)[,1];
qfD=0;

if(quadclassunit(D)[1]!=1,qfD = Listquadform(D), qfD = Listquadform1(D));

polform = qfD[3];

smatrix = matrix(2,#pN,i,j,squareN(D,pN[j]^factor(N)[j,2])[i]);
svec = matrix(#polform, #pN, i, j, vector( 2, k, if(Mod(Vec(qfD[1][i])[1],pN[j])!=Mod(0,pN[j]), [Mod((-Vec(qfD[1][i])[2]+smatrix[k,j])/(2*Vec(qfD[1][i])[1]),pN[j]^factor(N)[j,2]), Mod(1,pN[j]^factor(N)[j,2])], Polsol1(polform[i],pN[j]^factor(N)[j,2])[k]))); 

mats = matrix(#polform, 2^#pN,i,j, ChinSol1(vector(#pN,k,svec[i,k][binaryvector(#pN,j-1)[k]+1])));

matrep = matrix(#polform, 2^#pN, i, j, mattranspose(mats[i,j])*qfD[2][i]*mats[i,j]);

taurep = matrix(#polform, 2^#pN, i, j, (-matrep[i,j][1,2]*2+I*sqrt(-D))/(2*matrep[i,j][1,1]));

qtaurep = matrix(#polform, 2^#pN, i, j,exp(2*taurep[i,j]*I*Pi));

zrep = matrix(#polform, 2^#pN, i, j, subf1(qtaurep[i,j]));

Prep = matrix(#polform, 2^#pN, i, j, ellztopoint(E,zrep[i,j]));

zbyrep = vector(2^#pN);
Pbyrep = Vec(0,2^#pN);
for(i=1,2^#pN, for(j=1,#polform, zbyrep[i]+=zrep[j,i]));
for(i=1,2^#pN, Pbyrep[i] = ellztopoint(E,zbyrep[i]));


return([smatrix,svec,mats,matrep,taurep,qtaurep,zrep,Prep,zbyrep,Pbyrep]);

}


listrep(E,D,j)={

lists = List();
mats = Vec(0,2^#pN);
matrep = Vec(0,2^#pN);
taurep = Vec(0,2^#pN);
qtaurep = Vec(0,2^#pN);
zrep = Vec(0,2^#pN);
Prep = Vec(0,2^#pN);
zsum = 0;
Psum = 0;



for(i=0,2^#pN-1,listput(lists,List(sRep(polForm[j],binaryvector(#pN,i)))));

for(i=1,2^#pN,mats[i]=ChinSol(lists[i]));

for(i=1,2^#pN,matrep[i]=mattranspose(mats[i])*qfD[2][j]*mats[i]);

for(i=1,2^#pN, taurep[i]=(-matrep[i][1,2]*2+I*sqrt(-D))/(2*matrep[i][1,1]));

for(i=1,2^#pN, qtaurep[i] = exp(2*taurep[i]*I*Pi));

for(i=1,2^#pN, zrep[i] = subf1(qtaurep[i]));

for(i=1,2^#pN, zsum = zsum + zrep[i]);

for(i=1,2^#pN, Prep[i] = ellztopoint(E,zrep[i]));

Psum = ellztopoint(E,zsum);


return([matrep,taurep,qtaurep,zrep,Prep,zsum,Psum,mats]);


}


sumbyRep(E,vecrepr)={

zvecbyrep = Vec(0,2^#pN);
Pvecbyrep = Vec(0,2^#pN);
for(i=1,2^#pN, for(j=1,#polForm,zvecbyrep[i]+=vecrepr[j][4][i]));
for(i=1,2^#pN, Pvecbyrep[i] = ellztopoint(E,zvecbyrep[i]));

return([zvecbyrep,Pvecbyrep]);

}

CMpoints(L,quadformvec)={

return(0);

}

ChinSol(L)={

v = [Mod(1,1),Mod(1,1)];

for(i = 1, #L, v = [chinese(Mod(L[i][1],pN[i]^factor(N)[i,2]),v[1]), chinese(Mod(L[i][2],pN[i]^factor(N)[i,2]),v[2])]);

if(min(lift(v[1]),abs(lift(v[1])-N))==lift(v[1]), v11 = lift(v[1]), v11 = lift(v[1])-N);
v21 = lift(v[2]);

v12 = -bezout(v11,v21)[2];
v22 = bezout(v11,v21)[1];

M = [v11,v12;v21,v22];

return(M);

}

ChinSol1(L)={

v = [Mod(1,1),Mod(1,1)];

for(i = 1, #L, v = [chinese(L[i][1],v[1]), chinese(L[i][2],v[2])]);

if(min(lift(v[1]),abs(lift(v[1])-N))==lift(v[1]), v11 = lift(v[1]), v11 = lift(v[1])-N);
v21 = lift(v[2]);


v12 = -bezout(v11,v21)[2];
v22 = bezout(v11,v21)[1];

M = [v11,v12;v21,v22];

return(M);

}

squareN(D,N)={

s=[0,0];

for(i=0,N-1,if(Mod(i^2,N)==Mod(D,N),s[2]=i));

s[2]=s[2]-N;
s[1]=-s[2];

return(s);

}


orderPoint(E,vectorz)={

Pgen = ellgenerators(E)[1];
zgen = ellpointtoz(E,Pgen);

return(vector(#vectorz,z,if(!lindep([zgen,vectorz[z],E.omega[1],E.omega[2]])[2], lindep([zgen,vectorz[z],E.omega[1],E.omega[2]])[1]/lindep([zgen,vectorz[z],E.omega[1],E.omega[2]])[2]),lindep([zgen,vectorz[z],E.omega[1],E.omega[2]])[1]);

}

VectorDiscrim(n)={

VD = vector(2*n);
for(i=1,n,VD[2*i-1]=-4*i+1);
for(i=1,n,VD[2*i]=-4*i);

VDl = List(VD);

return(VD);

}

DiscHH(E,n)={

VD = VectorDiscrim(n);

DiscList = List();

for(i=1, #VD, if(HeegHy(E,VD[i])==0,listput(DiscList,VD[i])));

return(DiscList);

}

b0Heeg(E,n)={

Vd = DiscHH(E,n);

Pvector = List();

for(i=1 ,#Vd , listput(Pvector,orderPoint(E,sumbyRep(E,HeegnerPoints(E,Vd[i]))[1])));

return(Pvector);

}

b0Heeg1(E,n)={

Vd = DiscHH(E,n);

Pvector = List();

for(i=1 ,#Vd , listput(Pvector,orderPoint(E,HeegnerPoints2(E,Vd[i])[9])));

return(Pvector);

}




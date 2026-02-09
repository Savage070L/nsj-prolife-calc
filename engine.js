(function(){
var _z="eyJtIjpbMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC4wLDEuMzksMS43NSwxLjg0LDEuODQsMS44NSwxLjg3LDEuODgsMS44OCwxLjg5LDEuOTEsMS45MiwxLjkzLDEuOTMsMS45NCwxLjk2LDEuOTcsMS45NywxLjk4LDIuMCwyLjExLDIuMjcsMi40NCwyLjY2LDIuODksMy4xNCwzLjM5LDMuNjYsMy45Nyw0LjMsNC42OSw1LjE0LDUuNzIsNi4zNCw3LjAzLDcuNzQsOC40Niw5LjIxLDEwLjAsMTAuODgsMTEuODgsMTMuMDMsMTQuMzUsMTUuODEsMTcuNDEsMTkuMDcsMjAuODMsMjIuNzIsMjQuOCwyNy4xNCwyOS44LDMyLjgyLDM2LjIyLDQwLjA1LDQ0LjMzLDQ5LjExLDU0LjI3LDU5Ljk3LDY2LjI3LDczLjIzLDgwLjkyLDg5LjQyLDk4LjgyLDEwOS4yLDEyMC42NywxMzMuMzUsMTQ3LjM2LDE2Mi44NCwxNzkuOTUsMTk4Ljg2LDIxOS43NiwyNDIuODUsMjY4LjM3LDI5Ni41NywzMjcuNzMsMzYyLjE3LDQwMC4yMiw0NDIuMjcsNDg4Ljc0LDU0MC4wOSw1OTYuODQsNjU5LjU1LDcyOC44NSw4MDUuNDMsODkwLjA2LDk4My41OF0sImYiOlswLjAsMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMS4wNiwxLjA4LDAuOTgsMS4wMSwxLjAyLDEuMDQsMS4wNiwxLjA4LDEuMSwxLjExLDEuMTMsMS4xNSwxLjE3LDEuMiwxLjI0LDEuMjgsMS4zMywxLjQxLDEuNDgsMS41NSwxLjY2LDEuNzUsMS44NSwxLjk2LDIuMDYsMi4xOSwyLjM0LDIuNTEsMi43LDIuOTEsMy4xNCwzLjQsMy42OCw0LjAsNC4zMyw0LjY5LDUuMDYsNS40NCw1Ljg2LDYuMjksNi44LDcuNDYsOC4yMSw5LjA0LDkuOTMsMTAuODksMTEuOSwxMy4wNCwxNC4zMiwxNS43OSwxNy41LDE5LjQ1LDIxLjc2LDI0LjQ2LDI3LjY2LDMwLjk0LDM0LjYxLDM4LjcyLDQzLjMyLDQ4LjQ2LDU0LjIxLDYwLjY0LDY3Ljg0LDc1Ljg5LDg0LjksOTQuOTgsMTA2LjI1LDExOC44NiwxMzIuOTcsMTQ4Ljc1LDE2Ni40LDE4Ni4xNSwyMDguMjQsMjMyLjk1LDI2MC42LDI5MS41MywzMjYuMTMsMzY0LjgzLDQwOC4xMyw0NTYuNTcsNTEwLjc2LDU3MS4zOCw2MzkuMTksNzE1LjA1LDc5OS45MV0sImNtIjpbMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC4wLDAuMzIsMC4zMiwwLjMyLDAuMzIsMC4zMiwwLjM1LDAuNDMsMC40MywwLjQ1LDAuNDUsMC40NSwwLjQyLDAuNDUsMC40NSwwLjUzLDAuNjYsMC43NSwwLjkxLDEuMTIsMS4yOSwxLjU1LDEuODYsMi4yOCwyLjc3LDMuMzIsMy43NSw0LjM1LDQuNTgsNC44Myw1LjI3LDUuNzEsNi4xNCw2LjksNy43Miw4LjQyLDkuMzEsMTAuMzEsMTEuNTcsMTIuNTksMTMuNywxNS4wNCwxNi41OCwxNy42MSwxOC41NiwxOS4zOSwyMy4xNywyNy4wOSwyOC44MSwzMC41NywzMi40MSwzNS45NywzOS45Miw0NC4zLDQ5LjE2LDU0LjU2LDYwLjU1LDY3LjIsNzQuNTgsODIuNzcsOTEuODYsMTAxLjk0LDExMy4xMywxMjUuNTUsMTM5LjMzLDE1NC42MywxNzEuNjEsMTkwLjQ1LDIxMS4zNiwyMzQuNTYsMjYwLjMxLDI4OC44OSwzMjAuNiwzNTUuOCwzOTQuODYsNDM4LjIxLDQ4Ni4zMiw1MzkuNzEsNTk4Ljk2LDY2NC43MSw3MzcuNjgsODE4LjY2LDkwOC41MywxMDA4LjI3LDExMTguOTYsMTI0MS44XSwiY2YiOlswLjAsMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC4wLDAuMCwwLjAsMC40NiwwLjQ2LDAuNDYsMC40NiwwLjQ2LDAuNDgsMC42MiwwLjY4LDAuNzQsMC44MywwLjkzLDAuOTksMS4xMywxLjIyLDEuMywxLjQxLDEuNjIsMS43MywxLjkxLDIuMSwyLjMzLDIuNTUsMi43NywyLjk3LDMuMjIsMy41LDMuNzQsMy44Niw0LjEsNC40Miw0Ljc3LDUuMTIsNS41Nyw1Ljk5LDYuMzUsNi42LDcuMTksNy42OSw3Ljk3LDguMzIsOC44MSw5LjI3LDkuNjQsMTAuMzksMTAuOTksMTIuODMsMTQuNzgsMTUuNzksMTYuNiwxNy41LDE5LjIzLDIxLjEzLDIzLjIxLDI1LjUsMjguMDIsMzAuNzgsMzMuODIsMzcuMTYsNDAuODMsNDQuODYsNDkuMjksNTQuMTUsNTkuNDksNjUuMzYsNzEuODEsNzguODksODYuNjcsOTUuMjIsMTA0LjYxLDExNC45MywxMjYuMjcsMTM4LjczLDE1Mi40MiwxNjcuNDYsMTgzLjk4LDIwMi4xMywyMjIuMDcsMjQzLjk4LDI2OC4wNSwyOTQuNDksMzIzLjU0LDM1NS40NiwzOTAuNTMsNDI5LjA2LDQ3MS4zOV19";
var _q0=JSON.parse(atob(_z));
var _h0=0.07,_h1=0.12,_w0=100,_w2=0.01;
var _fa={annual:1.0,semiannual:0.51,quarterly:0.2575,monthly:0.0875};
var _ex={"3":{K:0.11923,L:0.02,M:0.021,N:0.02,O:0.0},"4":{K:0.1523,L:0.02,M:0.025,N:0.02,O:0.0},"5":{K:0.18538,L:0.04,M:0.027,N:0.025,O:0.0},"6":{K:0.2946,L:0.1434,M:0.029,N:0.03,O:0.003},"7":{K:0.3387,L:0.1623,M:0.030,N:0.03,O:0.003},"8":{K:0.3828,L:0.1812,M:0.03,N:0.03,O:0.003},"9":{K:0.4269,L:0.2001,M:0.03,N:0.03,O:0.003},"10":{K:0.501,L:0.249,M:0.03,N:0.04,O:0.003},"11":{K:0.5451,L:0.2679,M:0.03,N:0.04,O:0.003},"12":{K:0.5892,L:0.2868,M:0.03,N:0.04,O:0.003},"13":{K:0.6333,L:0.3057,M:0.03,N:0.04,O:0.003},"14":{K:0.6774,L:0.3246,M:0.03,N:0.04,O:0.003},"15":{K:0.7215,L:0.3435,M:0.025,N:0.04,O:0.003}};
var _rd={a:{t:0.001,e:0.25,q:0.4},b:{t:0.001,f:0.333333,e:0.25,q:0.4},c:{t:0.00047,e:0.25,q:0.4},d:{t:0.000141,e:0.25,q:0.4},e:{t:0.0045,e:0.05,q:0.25},f:{t:0.002,e:0.05,q:0.25},g:{t:0.00075,e:0.05,q:0.25}};

function _b1(r,i){
var v=1/(1+i),u=[],l=[],D=[],N=[],C=[],M=[];
var x;
for(x=0;x<101;x++)u[x]=r[x]/1000;
u[100]=1.0;
l[0]=1000000;
for(x=0;x<100;x++)l[x+1]=l[x]*(1-u[x]);
for(x=0;x<101;x++)D[x]=l[x]*Math.pow(v,x);
N[100]=D[100];
for(x=99;x>=0;x--)N[x]=D[x]+N[x+1];
for(x=0;x<100;x++)C[x]=(l[x]-l[x+1])*Math.pow(v,x+1);
C[100]=(l[100]>=1)?1:0;
M[100]=C[100];
for(x=99;x>=0;x--)M[x]=C[x]+M[x+1];
return{D:D,N:N,C:C,M:M,l:l,q:u};
}

function _b2(rd,rc,i){
var v=1/(1+i),qd=[],qc=[],x;
for(x=0;x<101;x++){qd[x]=rd[x]/1000;qc[x]=rc[x]/1000;}
qd[100]=1.0;
for(x=0;x<101;x++){
var s=qd[x]+qc[x];
if(s>1){qd[x]=qd[x]/s;qc[x]=qc[x]/s;}
}
var l=[],D=[],N=[],C=[],M=[],AG=[],Nq=[];
l[0]=1000000;
for(x=0;x<100;x++)l[x+1]=l[x]*(1-qd[x]-qc[x]);
for(x=0;x<101;x++)D[x]=l[x]*Math.pow(v,x);
N[100]=D[100];
for(x=99;x>=0;x--)N[x]=D[x]+N[x+1];
for(x=0;x<100;x++)C[x]=(l[x]-l[x+1])*Math.pow(v,x+1);
C[100]=(l[100]>=1)?1:0;
M[100]=C[100];
for(x=99;x>=0;x--)M[x]=C[x]+M[x+1];
for(x=0;x<101;x++)AG[x]=D[x]*qc[x];
Nq[100]=AG[100];
for(x=99;x>=0;x--)Nq[x]=AG[x]+Nq[x+1];
return{D:D,N:N,C:C,M:M,l:l,Nq:Nq};
}

function _ek(n){return String(Math.min(Math.max(n,3),15));}

function _gx(n,t){
var g6=_ex[_ek(n)].N,g7=_ex[_ek(n)].O,g10=Math.ceil(_ex[_ek(n)].M*100)/100,g2,g3;
if(t===1){g2=g10;g3=g10;}
else{var kt=_ek(t);g2=Math.ceil(_ex[kt].K*100)/100;g3=Math.ceil(_ex[kt].L*100)/100;}
return{g2:g2,g3:g3,g6:g6,g7:g7};
}

function _p1(tb,x,n,t){
var Ax=(tb.M[x]-tb.M[x+n]+tb.D[x+n])/tb.D[x];
var an=(tb.N[x]-tb.N[x+n])/tb.D[x];
var at=(tb.N[x]-tb.N[x+t])/tb.D[x];
var ge=_gx(n,t);
var BP=(Ax+ge.g7*an)/(at-ge.g6*at-(ge.g2+ge.g3*tb.D[x+1]/tb.D[x]));
return{Ax:Ax,an:an,at:at,BP:BP,ge:ge};
}

function _r1(tb,x,n,t,BP,ge,SA){
var rv=[],k,Ak,ak,atk,al,rr,sr,res,sur;
for(k=1;k<=n;k++){
Ak=(tb.M[x+k]-tb.M[x+n]+tb.D[x+n])/tb.D[x+k];
ak=(tb.N[x+k]-tb.N[x+n])/tb.D[x+k];
atk=(k<t)?(tb.N[x+k]-tb.N[x+t])/tb.D[x+k]:0;
al=(k===1)?ge.g3:0;
rr=Ak+ge.g7*ak-BP*(atk-ge.g6*atk-al);
sr=rr-(1-rr)*_w2;
res=rr*SA;
sur=Math.max(sr*SA,0);
if(k===n)sur=SA;
rv.push({year:k,reserve:Math.round(res),surrender:Math.round(sur)});
}
return rv;
}

function _p2a(rk,rs,n,ff,sg){
var rc=_rd[rk];
var bt=(rk==='b')?(rc.t*rc.f):rc.t;
var gt=Math.round((bt*(1+rc.e))/(1-rc.q)*10000)/10000;
var an=gt*rs;
if(sg)return{annual:Math.round(an),premium:Math.round(an*n),gross_tariff:gt};
return{annual:Math.round(an),premium:Math.round(an*ff),gross_tariff:gt};
}

function _p3(ctb,x,n,t,ge,ff,cs,sg){
var Ax=(ctb.M[x]-ctb.M[x+n]+ctb.D[x+n])/ctb.D[x]+(ctb.Nq[x]-ctb.Nq[x+n])/ctb.D[x];
var an=(ctb.N[x]-ctb.N[x+n])/ctb.D[x];
var at=(ctb.N[x]-ctb.N[x+t])/ctb.D[x];
var BP=(Ax+ge.g7*an)/(at-ge.g6*at-(ge.g2+ge.g3*ctb.D[x+1]/ctb.D[x]));
BP=Math.round(BP*10000)/10000;
var pm;
if(sg)pm=Math.round(BP*cs*n);
else pm=Math.round(BP*cs*ff);
return{BP:BP,premium:pm};
}

function _p4(tb,x,n,t,i,ff,ap,sg){
var v=1/(1+i),sm=0,k;
for(k=1;k<t;k++){
var ar=(tb.N[x+k]-tb.N[x+n])/tb.D[x+k];
sm+=ar*Math.pow(v,k);
}
var av=(t>1)?sm/(t-1):0;
var gt=(_rd.d.t*(1+0.25))/0.6;
var wa=av*gt*ap;
if(sg)return Math.round(wa*n);
return Math.round(wa*ff);
}

function _ca(p){
var db=new Date(p.dob),td=new Date(),ag=td.getFullYear()-db.getFullYear();
var m=td.getMonth()-db.getMonth();
if(m<0||(m===0&&td.getDate()<db.getDate()))ag--;
var x=ag,n=parseInt(p.term),fr=p.frequency||'annual',md=p.mode||'sa_to_premium';
var sg=(fr==='single'),i=(sg&&n===3)?_h1:_h0;
var gn=(p.gender==='male')?'m':'f';
var qx=_q0[gn],ff=sg?1:_fa[fr];
var tb=_b1(qx,i);
var t=sg?1:n;
var r=_p1(tb,x,n,t);
var BP=r.BP,ge=r.ge;
var SA,ap,gp;
if(md==='premium_to_sa'){
var pm=parseFloat(p.premium);
if(sg){SA=Math.round(pm/BP);}
else{SA=Math.round(pm/(BP*ff));}
ap=Math.round(BP*SA);
gp=sg?ap:Math.round(BP*SA*ff);
}else{
SA=parseFloat(p.sum_assured);
ap=Math.round(BP*SA);
gp=sg?ap:Math.round(BP*SA*ff);
}
var rv=_r1(tb,x,n,t,BP,ge,SA);
var ri={},rt=0;
if(p.riders){
var rds=p.riders;
if(rds.accidental_death&&rds.accidental_death.sum>0){
var o=_p2a('a',rds.accidental_death.sum,n,ff,sg);
ri.accidental_death=o;rt+=o.premium;
}
if(rds.traffic_death&&rds.traffic_death.sum>0){
var o=_p2a('b',rds.traffic_death.sum,n,ff,sg);
ri.traffic_death=o;rt+=o.premium;
}
if(rds.disability_accident&&rds.disability_accident.sum>0){
var o=_p2a('c',rds.disability_accident.sum,n,ff,sg);
ri.disability_accident=o;rt+=o.premium;
}
if(rds.disability_any&&rds.disability_any.sum>0){
var o=_p2a('d',rds.disability_any.sum,n,ff,sg);
ri.disability_any=o;rt+=o.premium;
}
if(rds.trauma&&rds.trauma.sum>0){
var o=_p2a('e',rds.trauma.sum,n,ff,sg);
ri.trauma=o;rt+=o.premium;
}
if(rds.temporary_disability&&rds.temporary_disability.sum>0){
var o=_p2a('f',rds.temporary_disability.sum,n,ff,sg);
ri.temporary_disability=o;rt+=o.premium;
}
if(rds.hospitalization&&rds.hospitalization.sum>0){
var o=_p2a('g',rds.hospitalization.sum,n,ff,sg);
ri.hospitalization=o;rt+=o.premium;
}
if(rds.critical_illness&&rds.critical_illness.sum>0){
var gci=(p.gender==='male')?'cm':'cf';
var ctb=_b2(qx,_q0[gci],i);
var o=_p3(ctb,x,n,t,ge,ff,rds.critical_illness.sum,sg);
ri.critical_illness=o;rt+=o.premium;
}
if(rds.premium_waiver&&rds.premium_waiver.enabled){
var wp=_p4(tb,x,n,t,i,ff,ap,sg);
ri.premium_waiver={premium:wp};rt+=wp;
}
}
return{
age:x,
term:n,
frequency:fr,
mode:md,
sum_assured:SA,
BP:Math.round(BP*1000000)/1000000,
annual_premium:ap,
gross_premium:gp,
reserves:rv,
riders:ri,
riders_total:rt,
total_premium:gp+rt
};
}

window.NSJEngine={calculate:_ca};
})();

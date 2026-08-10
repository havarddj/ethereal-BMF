
AttachSpec("spec");

d := 7;

K := QuadFld(d);
ZK := MaximalOrder(K);


level1 := LMFDBIdeal(K,"25.1");
level2 := LMFDBIdeal(K,"3025.1");

p := 2;

wt0 := BianchiWeight(K, 0,0);
wtp := BianchiWeight(K, 0,0 : char := p);

printf "Computing Bianchi cohomology space of level %o in characteristic %o\n", LMFDBLabel(level1),p;
B1 := BianchiCohomologySpace(level1,wtp);

printf "Computing Bianchi cohomology space of level %o in characteristic %o\n", LMFDBLabel(level2),p;
B2 := BianchiCohomologySpace(level2,wtp);

printf "Computing Bianchi cohomology space of level %o in characteristic 0\n", LMFDBLabel(level2);
B3 := BianchiCohomologySpace(level2,wt0);

E := Eigenforms(B1);

for e in E do 
	if e`eigenvalues[ideal<ZK | [1,2]>] eq 1 then
		f := e;
		break e;
	end if;
end for;

printf "Ethereal eigenvalue system in mod %o cohomology realised by %o\n", p, f;

D := level2/level1;
level_raise_f := [RaiseCocycleLevel(B1,B2,f,dd) : dd in Divisors(D)];


FF := sub<B2`forms | [u`vector : u in level_raise_f]>;
levelRaise := FF meet CharacteristicZeroImage(B2,B3);

printf "Dimension of level raise space: %o\n", Dimension(FF);
printf "Dimension of intersection of level raise space and characteristic 0 image: %o\n", Dimension(levelRaise); 

if Dimension(levelRaise) ge 1 then 
	for i in [1..Dimension(levelRaise)] do
		v := levelRaise.i;
		printf "Lift to characteristic zero realised by %o\n", Solution(Matrix([u`vector : u in level_raise_f]),v);
	end for;
end if;

P := LMFDBIdeal(K,"11.1");
aP := Eigenvalue(f,P);
printf "Eigenvalue of f at ideal with label %o: %o", LMFDBLabel(P), aP;
printf "a_P^2 - (1+N(P))^2 = %o\n", (aP^2 - (1+Norm(P))^2);
printf "P is a divisor of level*p? %o\n", level1*p subset P;


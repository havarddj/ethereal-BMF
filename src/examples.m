
load "look_for_lifts.m";
load "homology_torsion.m";


// Verifies Example 1.7
procedure Example_1_7()
	
	print "Verifying Example 1.7: eigenvalues over F_9 lifting to characteristic 0";
	str := "293.1;3;[ 0, 2 ];[ 2, 0 ];[ 0, 0 ];[ 2, 1 ];[ 1, 2 ];[ 1, 2 ];[ 1, 1 ];[ 1, 2 ];[ 2, 2 ];[ 0, 1 ];[ 0, 2 ];[ 2, 2 ];[ 2, 0 ];[ 1, 2 ];[ 1, 1 ];[ 2, 2 ];[ 1, 0 ];[ 2, 0 ];[ 0, 1 ];[ 0, 0 ];[ 0, 0 ];[ 1, 1 ];[ 0, 1 ];[ 2, 1 ];[ 0, 1 ];1;1;x^2 + 2*x + 2";
	f := LoadForm_nonEis(1,str);
	B := Parent(f);
	p := Characteristic(B);
	level := Level(B);
	FF<w> := CoefficientRing(f);
	SetPowerPrinting(FF,false);
	K := B`field;
	ZK := MaximalOrder(K);

	liftLevel := level*LMFDBIdeal(K,"5.2");
	printf "Computing mod %o space at level %o\n", p, LMFDBLabel(liftLevel);
	Bp := BianchiCohomologySpace(liftLevel, Weight(B));
	Bp := ChangeRing(Bp,FF);
	printf "Computing char. 0 space at level %o\n", LMFDBLabel(liftLevel);
	B0 := BianchiCohomologySpace(liftLevel,BianchiWeight(K,0,0));

	printf "Lifting f to higher level with degeneracy maps\n";
	lifts := sub<Bp`forms | [RaiseCocycleLevel(B,Bp,f,dd)`vector : dd in Divisors(liftLevel/level)]>;
	red := ChangeRing(CharacteristicZeroImage(Bp,B0),FF);

	assert Dimension(lifts meet red) gt 0;
	print "Verified existence of form in characteristic 0 lifting eigenvalues of form in Example 1.7";

	H := HeckeOperator(B0,LMFDBIdeal(K,"2.1"));
	F := NumberField(Polynomial([-1,1,1]));
	r1 := Roots(CharacteristicPolynomial(H));
	r2 := Roots(ChangeRing(CharacteristicPolynomial(H),F));
	printf "Hecke operator at 2.1 splits over Q: %o\n", &+[Integers()| u[2] : u in r1] eq Dimension(B0);
	printf "Hecke operator at 2.1 splits over Q(sqrt(5)): %o\n", &+[Integers()| u[2] : u in r2] eq Dimension(B0);
end procedure;


// Verifies Example 1.9
procedure Example_1_9()

	print "Verifying Example 1.9: Hecke field of lift is higher degree than that of ethereal form";
	str := "607.1;31;16;9;4;13;2;17;3;22;22;14;8;21;19;29;20;12;8;20;6;15;25;26;21;7;14;1;1;x + 30";
	f := LoadForm_nonEis(3,str);
	B := Parent(f);
	p := Characteristic(B);
	level := Level(B);
	K := B`field;
	ZK := MaximalOrder(K);

	liftLevel := level*LMFDBIdeal(K,"13.2");
	printf "Computing mod %o space at level %o\n", p, LMFDBLabel(liftLevel);
	Bp := BianchiCohomologySpace(liftLevel, Weight(B));
	printf "Computing char. 0 space at level %o\n", LMFDBLabel(liftLevel);
	B0 := BianchiCohomologySpace(liftLevel,BianchiWeight(K,0,0));

	printf "Lifting f to higher level with degeneracy maps\n";
	lifts := sub<Bp`forms | [RaiseCocycleLevel(B,Bp,f,dd)`vector : dd in Divisors(liftLevel/level)]>;
	red := CharacteristicZeroImage(Bp,B0);

	assert Dimension(lifts meet red) gt 0;
	print "Verified existence of form in characteristic 0 lifting eigenvalues of form in Example 1.9";

	H := HeckeOperator(B0,LMFDBIdeal(K,"3.1"));
	F := NumberField(Polynomial([-2,0,1]));
	r1 := Roots(CharacteristicPolynomial(H));
	r2 := Roots(ChangeRing(CharacteristicPolynomial(H),F));
	printf "Hecke operator at 3.1 splits over Q: %o\n", &+[Integers()| u[2] : u in r1] eq Dimension(B0);
	printf "Hecke operator at 3.1 splits over Q(sqrt(2)): %o\n", &+[Integers()| u[2] : u in r2] eq Dimension(B0);
end procedure;


// Verifies Example 2.5
procedure Example_2_5()

	print "Verifying Example 2.5: a large ethereal characteristic";
	K := QuadFld(7);
	ZK := MaximalOrder(K);
	level := LMFDBIdeal(K,"547.2");

	print "Computing torsion in H_1(Gamma0(547.2),Z)";
	AQinvs := ExtraModPClasses(level);

	printf "Torsion classes occur for primes in the set %o\n", Set(PrimeFactors(AQinvs[#AQinvs]));
end procedure;


// Verifies Example 3.9
procedure Example_3_9()

	print "Verifying Example 3.9: lift to an elliptic curve";
	str := "97.1;5;2;2;3;4;1;1;0;3;2;1;3;3;1;0;1;0;2;3;4;1;2;3;0;1;0;1;1;x + 4";
	f := LoadForm_nonEis(1,str);
	K := Parent(f)`field;
	ZK := MaximalOrder(K);

	E := EllipticCurve([1, K![-1,1], K![0,1], K![-4,-4], K![0,-4]]);
	PP := PrimesUpTo(100,K);

	printf "Eigenform of level %o:\n%o\n", LMFDBLabel(Level(f)), f;
	printf "Elliptic curve with LMFDB label 2.0.4.1-194.1-a1:\n%o\n", E;

	printf "Traces of Frobenius match Hecke eigenvalues up to norm 100: %o\n", &and[Eigenvalue(f,P) eq TraceOfFrobenius(E,P) : P in PP | GCD(P,Conductor(E)) eq 1*ZK];
end procedure;


// Verifies Example 3.10
procedure Example_3_10()

	print "Verifying Example 3.10: probable lift not verifiable with cohomology";
	str := "289.1;2;0;0;0;0;1;1;1;1;1;0;0;0;0;1;1;0;0;0;0;0;0;0;0;0;0;1;1;x + 1";
	f := LoadForm_nonEis(3,str);
	B := Parent(f);
	p := Characteristic(B);
	level := Level(B);
	K := B`field;
	ZK := MaximalOrder(K);

	liftLevel := level*ideal<ZK | 2>^2;
	printf "Computing mod %o space at level %o\n", p, LMFDBLabel(liftLevel);
	Bp := BianchiCohomologySpace(liftLevel, Weight(B));
	printf "Computing char. 0 space at level %o\n", LMFDBLabel(liftLevel);
	B0 := BianchiCohomologySpace(liftLevel,BianchiWeight(K,0,0));

	printf "Lifting f to higher level with degeneracy maps\n";
	lifts := sub<Bp`forms | [RaiseCocycleLevel(B,Bp,f,dd)`vector : dd in Divisors(liftLevel/level)]>;
	red := CharacteristicZeroImage(Bp,B0);

	printf "Reduction intersects level-raise space non-trivially: \n", Dimension(lifts meet red) gt 0;
end procedure;


// Verifies Example 4.5
procedure Example_4_5()

	print "Verifying Example 4.5: ethereal eigenvalue system with no small level-raising primes";
	str := "421.1;17;4;10;7;3;11;8;6;5;11;6;7;15;16;6;13;12;4;5;10;1;11;1;8;1;1;x + 16";
	f := LoadForm_nonEis(7,str);
	B := Parent(f);
	p := Characteristic(B);
	level := Level(B);
	K := B`field;
	ZK := MaximalOrder(K);

	printf "Eigenform:\n%o\n", f;
	print "f satisfies the level-raising condition at prime P:";
	PP := [P : P in PrimesUpTo(100,K) | GCD(P,Level(B)*p) eq 1*ZK];
	print "Label | a_p | (a_p)^2 - (1+Norm(P))^2";
	for P in PP do 
		printf "%o | %o | %o\n", LMFDBLabel(P), Eigenvalue(f,P), Eigenvalue(f,P)^2 - (1+Norm(P))^2;
	end for;

	printf "Primes for which f satisfies the level-raising condition up to norm 100: \n%o\n", [LMFDBLabel(P) : P in PP | Eigenvalue(f,P)^2 - (1+Norm(P))^2 eq 0];
end procedure;


// Verifies Example 4.6
procedure Example_4_6()

	print "Verifying Example 4.6: lift with large Hecke field";
	str := "667.2;53;48;43;45;12;34;8;52;36;41;46;52;41;21;41;21;47;46;43;20;7;30;22;33;1;1;x + 52";
	f := LoadForm_nonEis(7,str);
	B := Parent(f);
	p := Characteristic(B);
	level := Level(B);
	K := B`field;
	ZK := MaximalOrder(K);

	printf "Ethereal form:\n%o\n", f;

	liftLevel := level*ideal<ZK | 2*K.1-1>;
	printf "Computing mod %o space at level %o\n", p, LMFDBLabel(liftLevel);
	Bp := BianchiCohomologySpace(liftLevel, Weight(B));
	printf "Computing char. 0 space at level %o\n", LMFDBLabel(liftLevel);
	B0 := BianchiCohomologySpace(liftLevel,BianchiWeight(K,0,0));

	printf "Lifting f to higher level with degeneracy maps\n";
	lifts := sub<Bp`forms | [RaiseCocycleLevel(B,Bp,f,dd)`vector : dd in Divisors(liftLevel/level)]>;
	red := CharacteristicZeroImage(Bp,B0);

	assert Dimension(lifts meet red) gt 0;
	print "Verified existence of form in characteristic 0 lifting eigenvalues of form in Example 4.6";

	P := LMFDBIdeal(K,"2.1");
	printf "Computing Hecke operator at prime %o\n", LMFDBLabel(P);
	H := HeckeOperator(B0,P);
	printf "Factorization of characteristic polynomial of Hecke operator at prime %o: \n%o\n", LMFDBLabel(P), Factorization(CharacteristicPolynomial(ChangeRing(H,Rationals())));

end procedure;


// Verifies Example 4.7
procedure Example_4_7()

	print "Verifying Example 4.7: lifting a mod 379 ethereal form to characteristic 0";
	str := "641.1;379;3;45;217;68;217;235;11;171;219;300;242;291;17;141;320;83;128;285;85;283;103;103;87;1;1;x + 378";
	f := LoadForm_nonEis(7,str);
	B := Parent(f);
	p := Characteristic(B);
	level := Level(B);
	K := B`field;
	ZK := MaximalOrder(K);

	liftLevel := level*LMFDBIdeal(K,"2.1");
	printf "Computing mod %o space at level %o\n", p, LMFDBLabel(liftLevel);
	Bp := BianchiCohomologySpace(liftLevel, Weight(B));
	printf "Computing char. 0 space at level %o\n", LMFDBLabel(liftLevel);
	B0 := BianchiCohomologySpace(liftLevel,BianchiWeight(K,0,0));

	printf "Lifting f to higher level with degeneracy maps\n";
	lifts := sub<Bp`forms | [RaiseCocycleLevel(B,Bp,f,dd)`vector : dd in Divisors(liftLevel/level)]>;
	red := CharacteristicZeroImage(Bp,B0);

	assert Dimension(lifts meet red) gt 0;
	print "Verified existence of form in characteristic 0 lifting eigenvalues of form in Example 4.7";
end procedure;


// Verifies Example 4.9
procedure Example_4_9()
	K := QuadFld(3);
	ZK := MaximalOrder(K);
	level := LMFDBIdeal(K,"931.2");
	p := 7;
	printf "Computing Bianchi cohomology space of level %o in characteristic %o\n", LMFDBLabel(level), p;
	B := BianchiCohomologySpace(level,BianchiWeight(K,0,0 : char := p));

	print "Computing Hecke operators up to norm 100";
	HP := GoodHeckePrimes(B,100);
	HH := [HeckeOperator(B,P) : P in HP];

	E := Eigenforms(B);

	printf "Number of eigenforms: %o\n", #E;
	printf "Dimensions of each form's eigenspace: %o\n", [Dimension(Eigenspace(ee)) : ee in E];
end procedure;


all_examples := [
	Example_1_7,
	Example_1_9,
	Example_2_5,
	Example_3_9,
	Example_3_10,
	Example_4_5,
	Example_4_6,
	Example_4_7,
	Example_4_8
];

for eg in all_examples do 
	eg();
	print "";
end for;

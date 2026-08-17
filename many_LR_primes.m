
load "look_for_lifts.m";

// d gives the imaginary quadratic field
// line is the line number of the form in the file d*_unliftable.csv
procedure FindLRPrimes(d,line : bound := 1000);

	rr := Read("sage/d" cat Sprint(d) cat "_unliftable.csv");
	ll := Split(rr,"\n");
	str := ll[line];
	f := LoadForm_liftable(d,str[1..#str-1] cat ";[]");
	B := Parent(f);
	K := B`field;
	ZK := B`integers;

	filename := "data/LR/" cat Sprint(d) cat "/" cat Split(str,";")[1] cat "_" cat Split(str,";")[2];

	LR_primes := [];
	for P in PrimesUpTo(bound,K) do 
		if GCD(P,Level(B)*Characteristic(B)) eq 1*ZK then
			e := Eigenvalue(f,P);
			if e^2 - (1+Norm(P)) eq 0 then 
				Append(~LR_primes,P);
				printf "Found a level-raise prime: %o\n", LMFDBLabel(P);
				Write(filename,LMFDBLabel(P));
			else 
				printf "Prime %o was not level-raise\n", LMFDBLabel(P);
			end if;
		end if;
	end for;

end procedure;



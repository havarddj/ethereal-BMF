
AttachSpec("../spec");

load "look_for_lifts.m";

d := 1;

K := QuadFld(d);
ZK := MaximalOrder(K);

rr := Read("sage/d" cat Sprint(d) cat "_liftable.csv");
lines := Split(rr);
filename := "verified_d" cat Sprint(d) cat ".csv";
SetColumns(0);

for str in lines[2..#lines] do 
	levelLabel, p, eigs := Explode(Split(str,";"));
	tt, liftLabel := VerifyLMFDBLift(d,str);
	if tt then 
		outstr := levelLabel cat ";" cat p cat ";" cat liftLabel cat ";" cat eigs;
		Write(filename,outstr);
		printf "Verified lift %o\n", str;
	else 
		printf "Failed to verify lift %o\n", str;
	end if;
end for;


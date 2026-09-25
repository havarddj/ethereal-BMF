load "look_for_lifts.m";

d := 7;

K := QuadFld(d);
ZK := MaximalOrder(K);

rr := Read("sage/d" cat Sprint(d) cat "_liftable.csv");
num := #Split(rr);

inds := [92];

for i in inds do 
	VerifyLift(d,i);
end for;


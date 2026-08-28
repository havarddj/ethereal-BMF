
// tells you whether e "blocks" lifting to an elliptic curve.
// i.e. whether there can be an integer n with abs(n) < 2*N(P)^(1/2)
IsBlockingEigenvalue := function(P,e)
	c := Floor(2*Norm(P)^(1/2));
	poss := [-c..c];
	return not 0 in [(e - u) : u in poss];
end function;

// tells us if e satisfies the level raise condition
IsLevelRaise := function(P,e) 
	return ((1+Norm(P))^2 - e^2) eq 0;
end function;


// if this is true for any of the eigenvalues in our range, we don't expect an EC lift 
BlockAndNotLR := function(P,e)
	return IsBlockingEigenvalue(P,e) and not IsLevelRaise(P,e);
end function;


load "look_for_lifts.m";


// prints data about blocking eigenvalues in lmfdbNonlift_d*.csv
// print_lift: whether to print those that should have a rational lift 
SurveyBlockingEigenvalues := function(d : noBlock := false)
	K := QuadFld(d);
	ZK := MaximalOrder(K);

	input := "data/lmfdbNonlift_d" cat Sprint(d) cat ".csv";
	lines := Split(Read(input), "\n");
	topLine := lines[1];
	topElts := Split(topLine, ";");
	primeLabels := topElts[3..#topElts -3];
	primes := [LMFDBIdeal(K,l) : l in primeLabels];

	// we count the number of rational systems that didnt have LMFDB lifts 
	rat_noLift := [];
	// and the number of these that are blocked by an eigenvalue a prime for which the system does not satisfy
	// the level-raising condition 
	rat_blocked := [];
	rat_blocked_lt20 := [];

	rat_not_blocked := [];
	rat_not_blocked_lt20 := [];


	for l in [2..#lines] do 

		ll := Split(lines[l],";");
		p := StringToInteger(ll[2]);
		F := GF(p);

		try
			eigs := [F!StringToInteger(u) : u in ll[3..#ll-3]];
			if #[i : i in [1..#eigs] | BlockAndNotLR(primes[i],eigs[i])] eq 0 then 
				if noBlock then 
					print lines[l];
				end if;
				Append(~rat_not_blocked,lines[l]);
				if p lt 20 then 
					Append(~rat_not_blocked_lt20,lines[l]);
				end if;
			else
				Append(~rat_blocked,lines[l]);
				if p lt 20 then 
					Append(~rat_blocked_lt20,lines[l]);
				end if;
			end if;
			Append(~rat_noLift,lines[l]);
		catch e 
			;
		end try;
	end for;


	return #rat_noLift - #rat_not_blocked_lt20, #rat_noLift - #rat_not_blocked;

end function;




Blah := function(d)

	K := QuadFld(d);
	ZK := MaximalOrder(K);


	input := "data/nonEis_d" cat Sprint(d) cat "_v2.csv";
	lines := Split(Read(input), "\n");
	topLine := lines[1];
	topElts := Split(topLine, ";");
	primeLabels := topElts[3..#topElts -3];
	primes := [LMFDBIdeal(K,l) : l in primeLabels];

	// all rational systems 
	eth_rat := [];

	// all rational systems with no blocking eigenvalues 
	eth_rat_noBlock := [];

	// same as above but for p at most 20 
	eth_rat_lt20 := [];
	eth_rat_lt20_noBlock := [];

	noBlock := false;

	for l in [2..#lines] do 
		ll := Split(lines[l],";");
		p := StringToInteger(ll[2]);
		F := GF(p);

		try
			eigs := [F!StringToInteger(u) : u in ll[3..#ll-3]];
			if #[i : i in [1..#eigs] | BlockAndNotLR(primes[i],eigs[i])] eq 0 then 
				if noBlock then 
					print lines[l];
				end if;
				Append(~eth_rat_noBlock,lines[l]);
				if p lt 20 and p gt 2 then
					Append(~eth_rat_lt20_noBlock,lines[l]);
				end if;
			end if;
			Append(~eth_rat,lines[l]);
			if p lt 20 and p gt 2 then 
				Append(~eth_rat_lt20,lines[l]);
			end if;
			
		catch e 
			;
		end try;
	end for;


	// now we load the unliftables, and check if the noBlocks are in there 
	input := "data/lmfdbNonlift_d" cat Sprint(d) cat ".csv";
	nonlift := Split(Read(input), "\n");
		
	// rational ethereal systems with no blocking eigenvalue that are LMFDB liftable 
	eth_rat_noBlock_lmfdb := [];
	// same but p < 20
	eth_rat_lt20_noBlock_lmfdb := [];
	
	for e in eth_rat_noBlock do 
		if not e in nonlift then 
			Append(~eth_rat_noBlock_lmfdb,e);
			p := StringToInteger(Split(e,";")[2]);
			if p lt 20 and p gt 2 then 
				Append(~eth_rat_lt20_noBlock_lmfdb,e);
			end if;
		end if;
	end for;
	

	//print eth_rat_lt20_noBlock;
	return #eth_rat_lt20_noBlock, #eth_rat_lt20_noBlock_lmfdb, #eth_rat_noBlock, #eth_rat_noBlock_lmfdb;
end function;


for d in [1,2,3,7,11] do 
	Blah(d);
end for;




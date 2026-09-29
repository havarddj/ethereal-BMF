

AttachSpec("BMF/spec");
load "loading.m";

// returns true if cohomology can verify the lift specified by the string `str`,
// a line in irrat_lifts_d*.csv.
function VerifyIrrationalLift(d,str)
	K := QuadFld(d);
	ZK := MaximalOrder(K);
	ss := Split(str,";");

	levelLabel := ss[1];
	level := LMFDBIdeal(K,levelLabel);
	liftLevel := LMFDBIdeal(K,ss[3]);
	p := ss[2];
	pp := StringToInteger(p);

	// we need to build the form from the data provided. that means find it in the file nonEis_d*.csv, 
	// by finding all forms with the right level and prime. if there is more than one of these, we 
	// need to check the eigenvalues match.

	nonEis := Split(Read("data/nonEis_d" cat Sprint(d) cat ".csv"));

	// these are the actual ideals represented by our LMFDB label strings 
	nn := #Split(nonEis[1],";");
	idealStrings := Split(nonEis[1],";")[3..nn-3];
	ideals := [LMFDBIdeal(K,u) : u in idealStrings];

	possible := [];
	for e in nonEis do 
		se := Split(e,";");
		if se[1] eq levelLabel and se[2] eq p then 
			Append(~possible,e);
		end if;
	end for;

	if #possible gt 1 then 
		// this case probably won't happen. but just in case! 
		_<x> := PolynomialRing(Rationals());
		f := eval ss[#ss];
		F := NumberField(f);
		ZF := MaximalOrder(F);
		char0Eigs := [F!(eval u) : u in ss[4..#ss-1]];

		// this collects all possible mod p reductions of the char 0 system 
		PP := [f[1] : f in Factorization(pp*ZF)];
		downs := [];
		for P in PP do 
			_,dd := ResidueClassField(P);
			Append(~downs,dd);
		end for;
		reds := [[d(u) : u in char0Eigs] : d in downs];

		// now we look for one that matches away from the lifted level and characteristic 
		matching := [];
		for l in possible do
			ll := Split(l,";");
			eigs := [StringToInteger(u) : u in ll[3..#ll-3]];

			for r in reds do 
				keep := true; 
				for i in [1..#ideals] do 
					if GCD(ideals[i],liftLevel*pp) eq 1*ZK then
						keep and:= (eigs[i]-r[i]) eq 0;
					end if;
				end for;
				if keep then
					Append(~matching,l);
				end if;
			end for;
		end for;

		if #matching gt 1 then 
			print "More than one form matching string:";
			print str;
			print "Please try again manually.";
			return false, "";
		else
			f := LoadForm_nonEis(d,matching[1]);
			form_str := matching[1];
		end if;
	else 
		f := LoadForm_nonEis(d,possible[1]);
		form_str := possible[1];
	end if;


	// now we do the usual level-raising procedure.
	B := Parent(f);

	// for annoying reasons we need to re-define K here. bah! 
	K := B`field;
	level := LMFDBIdeal(K,levelLabel);
	liftLevel := LMFDBIdeal(K,ss[3]);
	B0 := BianchiCohomologySpace(liftLevel, BianchiWeight(K,0,0));
	Bp := BianchiCohomologySpace(liftLevel, Weight(B));

	// when f is defined over a larger field than GF(p), we need to change all the rings.
	Rf := CoefficientRing(f);
	BBp := ChangeRing(Bp,Rf);

	levelRaise := [RaiseCocycleLevel(B,BBp,f,dd) : dd in Divisors(liftLevel/level)];
	red := ChangeRing(CharacteristicZeroImage(Bp,B0),Rf);
	LR := sub<BBp`forms | [u`vector : u in levelRaise]>;

	return Dimension(LR meet red) gt 0, form_str;
end function;



procedure VerifyIrrational(d)
	outFileName := "data/verified_irrational_d" cat Sprint(d) cat ".csv";
	lines := Split(Read("data/irrat_lifts_d" cat Sprint(d) cat ".csv"));
	SetColumns(0);
	for l in lines[2..#lines] do 
		ll := Split(l,";");
		if ll[3] ne "None found" and ll[3] ne "No level raising primes available" then
			tt, str := VerifyIrrationalLift(d,l);

			// if we have a verified lift, we build the string and add it to the relevant file
			if tt then 
				eigstr := "[";
				split := Split(str,";");
				for i in [3..#split-3] do
					eigstr cat:= split[i];
					if i ne #split-3 then
						eigstr cat:= ", ";
					else 
						eigstr cat:= "]";
					end if;
				end for;
				outStr := ll[1] cat ";" cat ll[2] cat ";" cat ll[3] cat ";" cat eigstr;
				print "Wrote verified lift to file!: " cat outStr;
				Write(outFileName,outStr);
			end if; 
		end if;
	end for;
end procedure;


